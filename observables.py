"""Safe, reusable optical observables for convergence, sweeps, and design."""

from __future__ import annotations

import ast
import math
import re
import statistics


_OBSERVABLE_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_]{0,31}$")
_RESERVED_NAMES = {"R", "T", "A", "R0", "T0", "Ro", "To", "abs", "sqrt",
                   "mean", "stdev", "min", "max"}


def _order_value(result: dict, port: str, m: int, n: int) -> float:
    if port not in ("R", "T") or not isinstance(m, int) or not isinstance(n, int):
        raise ValueError("Order functions require R or T and integer (m,n)")
    if abs(m) > 50 or abs(n) > 50:
        raise ValueError("Diffraction-order indices must be between -50 and 50")
    for row in result.get("orders", []):
        if int(row["m"]) == m and int(row["n"]) == n:
            return float(row[port])
    return 0.0


def evaluate_observable(result: dict, expression: str,
                        named_values: dict[str, float] | None = None) -> float:
    """Evaluate a restricted arithmetic expression against one scattering result.

    Scalars R, T, A, R0, T0 and order functions Ro(m,n), To(m,n) are
    available. Missing diffraction orders have zero far-field efficiency.
    Supported helpers are abs, sqrt, mean, stdev, min, and max.
    """
    expression = str(expression).strip()
    if not expression or len(expression) > 500:
        raise ValueError("Observable expression must contain 1–500 characters")
    try:
        tree = ast.parse(expression, mode="eval")
    except SyntaxError as exc:
        raise ValueError("Observable expression is not valid arithmetic") from exc
    if sum(1 for _ in ast.walk(tree)) > 120:
        raise ValueError("Observable expression is too complex")
    scalars = {key: float(result[key]) for key in ("R", "T", "A", "R0", "T0")}
    for key, value in (named_values or {}).items():
        if not _OBSERVABLE_NAME.fullmatch(str(key)) or key in _RESERVED_NAMES:
            raise ValueError(f"Invalid named observable: {key}")
        scalars[str(key)] = float(value)

    def number(node) -> float:
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            value = float(node.value)
        elif isinstance(node, ast.Name) and node.id in scalars:
            value = scalars[node.id]
        elif isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = number(node.operand) * (-1 if isinstance(node.op, ast.USub) else 1)
        elif isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)):
            left, right = number(node.left), number(node.right)
            if isinstance(node.op, ast.Add): value = left + right
            elif isinstance(node.op, ast.Sub): value = left - right
            elif isinstance(node.op, ast.Mult): value = left * right
            elif isinstance(node.op, ast.Div):
                if abs(right) < 1e-15: raise ValueError("Observable expression divides by zero")
                value = left / right
            else:
                if abs(right) > 8: raise ValueError("Observable exponents are limited to ±8")
                value = left ** right
        elif isinstance(node, (ast.List, ast.Tuple)):
            return [number(item) for item in node.elts]
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            name = node.func.id
            if node.keywords:
                raise ValueError("Observable functions do not accept named arguments")
            if name in ("Ro", "To"):
                if len(node.args) != 2:
                    raise ValueError(f"{name} needs integer m and n")
                indices = []
                for item in node.args:
                    value_index = number(item)
                    if value_index != int(value_index):
                        raise ValueError(f"{name} order indices must be integers")
                    indices.append(int(value_index))
                value = _order_value(result, name[0], *indices)
            elif name in ("abs", "sqrt") and len(node.args) == 1:
                arg = number(node.args[0])
                if name == "sqrt" and arg < 0:
                    raise ValueError("sqrt requires a nonnegative value")
                value = abs(arg) if name == "abs" else math.sqrt(arg)
            elif name in ("mean", "stdev", "min", "max"):
                args = [number(item) for item in node.args]
                values = args[0] if len(args) == 1 and isinstance(args[0], list) else args
                if not values:
                    raise ValueError(f"{name} needs at least one value")
                if name == "mean": value = statistics.fmean(values)
                elif name == "stdev": value = statistics.pstdev(values)
                elif name == "min": value = min(values)
                else: value = max(values)
            else:
                raise ValueError("Allowed functions: Ro, To, abs, sqrt, mean, stdev, min, max")
        else:
            raise ValueError("Observable expression contains an unsupported operation")
        if isinstance(value, list):
            return value
        if not math.isfinite(value):
            raise ValueError("Observable expression produced a non-finite value")
        return float(value)

    value = number(tree.body)
    if isinstance(value, list):
        raise ValueError("Observable expression must return one number")
    return value


def observable_spec(item: dict) -> dict:
    """Normalize a standard quantity or custom-expression specification."""
    quantity = str(item.get("quantity", "")).strip()
    if quantity in ("R", "T", "A", "R0", "T0"):
        return {"quantity": quantity, "expression": quantity,
                "label": str(item.get("label") or quantity)}
    expression = str(item.get("expression", "")).strip()
    if quantity not in ("expression", "custom") or not expression:
        raise ValueError("Choose R, T, A, R0, T0, or provide a custom observable expression")
    return {"quantity": "expression", "expression": expression,
            "label": str(item.get("label") or expression)[:120]}


def observable_value(result: dict, item: dict) -> float:
    return evaluate_observable(result, observable_spec(item)["expression"])


def normalize_observable_definitions(items: list[dict] | None) -> list[dict]:
    """Validate ordered reusable observable definitions.

    Later expressions may reference names defined earlier in the list.  The
    physical values are still evaluated independently for every wavelength.
    """
    if len(items or []) > 20:
        raise ValueError("Use no more than 20 named observables")
    normalized, seen = [], set()
    for item in items or []:
        name = str(item.get("name", "")).strip()
        expression = str(item.get("expression", "")).strip()
        if (not _OBSERVABLE_NAME.fullmatch(name) or name in _RESERVED_NAMES or
                name in seen):
            raise ValueError("Observable names must be unique identifiers and cannot replace built-in quantities")
        if not expression:
            raise ValueError(f"Named observable {name} needs an expression")
        # Parse now; unknown names are checked during ordered evaluation.
        try:
            ast.parse(expression, mode="eval")
        except SyntaxError as exc:
            raise ValueError(f"Observable {name} is not valid arithmetic") from exc
        expanded = expression
        for prior in normalized:
            expanded = re.sub(rf"\b{re.escape(prior['name'])}\b",
                              f"({prior['expanded_expression']})", expanded)
        normalized.append({"name": name, "expression": expression,
                           "expanded_expression": expanded,
                           "label": str(item.get("label") or name)[:120]})
        seen.add(name)
    return normalized


def evaluate_named_observables(result: dict, definitions: list[dict] | None) -> dict[str, float]:
    """Evaluate an ordered observable registry for one scattering result."""
    values: dict[str, float] = {}
    for item in normalize_observable_definitions(definitions):
        try:
            values[item["name"]] = evaluate_observable(result, item["expression"], values)
        except ValueError as exc:
            raise ValueError(f"Observable {item['name']}: {exc}") from exc
    return values


def design_metric_spec(item: dict, definitions: list[dict] | None = None) -> dict:
    """Resolve a built-in, inline, or named metric used by design tools."""
    registry = {row["name"]: row for row in normalize_observable_definitions(definitions)}
    named = str(item.get("observable", "")).strip()
    if named:
        if named not in registry:
            raise ValueError(f"Unknown named observable: {named}")
        row = registry[named]
        return {"quantity": "expression", "observable": named,
                "expression": row["expanded_expression"], "label": row["label"]}
    return observable_spec(item)


def design_metric_value(result: dict, item: dict, definitions: list[dict] | None = None,
                        named_values: dict[str, float] | None = None) -> float:
    """Evaluate one design metric while honoring a reusable registry."""
    spec = design_metric_spec(item, definitions)
    values = named_values if named_values is not None else evaluate_named_observables(result, definitions)
    if spec.get("observable"):
        return float(values[spec["observable"]])
    return evaluate_observable(result, spec["expression"], values)
