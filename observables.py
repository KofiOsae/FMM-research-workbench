"""Safe, reusable optical observables for convergence, sweeps, and design."""

from __future__ import annotations

import ast
import math
import statistics


def _order_value(result: dict, port: str, m: int, n: int) -> float:
    if port not in ("R", "T") or not isinstance(m, int) or not isinstance(n, int):
        raise ValueError("Order functions require R or T and integer (m,n)")
    if abs(m) > 50 or abs(n) > 50:
        raise ValueError("Diffraction-order indices must be between -50 and 50")
    for row in result.get("orders", []):
        if int(row["m"]) == m and int(row["n"]) == n:
            return float(row[port])
    return 0.0


def evaluate_observable(result: dict, expression: str) -> float:
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
