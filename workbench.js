/* Progressive enhancement; original operation handlers remain the single source of actions. */
(() => {
  const nav = document.querySelector('.operation-nav');
  const appearance=document.createElement('details');appearance.className='plot-appearance';appearance.innerHTML='<summary>Plot appearance · palette and smoothing</summary>';const appearanceControls=document.getElementById('vizSmoothing').closest('.row3');appearanceControls.before(appearance);appearance.append(appearanceControls);
  nav.insertAdjacentHTML('beforeend', '<div><label for="quickAction">Run from here</label><select id="quickAction"></select></div><button id="quickRun" type="button">Run selected</button><div><label for="resolutionLimit">My maximum scan points</label><input id="resolutionLimit" type="number" min="9" max="40401" value="10000"></div><a href="/FIRST_STEPS.html" target="_blank">First steps ↗</a>');
  const action = document.getElementById('quickAction'), run = document.getElementById('quickRun');
  function refreshActions() {
    const old = action.value;
    const buttons = [...document.querySelectorAll('main button[id]')].filter(b => !b.closest('[hidden]') && (/^run/i.test(b.id)||/^(calculate|find|fit|detect|adaptive|compute|validate)/i.test(b.textContent.trim())) && !/benchmark/i.test(b.id));
    action.replaceChildren(...buttons.map(b => {const o=document.createElement('option');o.value=b.id;o.textContent=b.textContent.trim();o.disabled=b.disabled;return o;}));
    if([...action.options].some(o=>o.value===old))action.value=old;
    run.disabled=!action.value || document.getElementById(action.value)?.disabled;
  }
  run.onclick=()=>document.getElementById(action.value)?.click();
  action.onchange=()=>{run.disabled=document.getElementById(action.value)?.disabled;};
  document.querySelector('.workspace-nav').addEventListener('click',()=>setTimeout(refreshActions,0));
  nav.addEventListener('change',()=>setTimeout(refreshActions,0));
  const observer=new MutationObserver(refreshActions);
  observer.observe(document.querySelector('main'),{subtree:true,attributes:true,attributeFilter:['hidden','disabled']});
  refreshActions();
  const stickySize=new ResizeObserver(()=>{const height=document.querySelector('.workspace-nav').getBoundingClientRect().height;nav.style.top=innerWidth>650?height+'px':'0';document.documentElement.style.setProperty('--study-top',(height+nav.getBoundingClientRect().height+12)+'px');});stickySize.observe(document.querySelector('.workspace-nav'));stickySize.observe(nav);
  const budgetNote=document.createElement('small');budgetNote.className='scan-estimate';nav.append(budgetNote);
  function estimateRun(){const id=action.value;let count=1,op='';const n=id=>Number(document.getElementById(id)?.value||1);if(id==='runScatteringMap'){count=n('mapLambdaPoints')*n('mapThetaPoints');op='angle_wavelength';}else if(id==='run'||id==='runTmmSpectrum'){count=n('points');op=id==='run'?'spectrum':'tmm';}else if(id==='runPolarizationMap'){count=n('polMapPoints')**2;op='polarization_kspace';}let seconds=null;try{seconds=Number(localStorage.getItem('timing:'+op))*count||null;}catch{}budgetNote.textContent=count>1?`${count.toLocaleString()} requested samples · ${seconds?'estimated '+Math.ceil(seconds)+' s from an earlier run':'time estimate after first completed samples'}`:'Run monitor shows elapsed time; solver stages without measurable work use an indeterminate bar.';}
  action.addEventListener('change',estimateRun);document.querySelector('main').addEventListener('input',estimateRun);estimateRun();

  const tray=document.createElement('div');tray.className='run-tray';tray.setAttribute('aria-live','polite');document.body.append(tray);
  const originalFetch=window.fetch.bind(window);
  const excluded=new Set(['materials','figure','import_material']);
  window.fetch=async function(input,options={}) {
    const url=String(input), operation=url.startsWith('/api/')?url.slice(5):'';
    if(options.method!=='POST'||!operation||excluded.has(operation)||operation.includes('/'))return originalFetch(input,options);
    let payload;try{payload=JSON.parse(options.body);}catch{return originalFetch(input,options);}
    if(operation==='vector_modes'&&document.getElementById('vmCoreShape')){payload.core_shape=document.getElementById('vmCoreShape').value;payload.sidewall_angle_deg=Number(document.getElementById('vmSidewall').value);options={...options,body:JSON.stringify(payload)};}
    const count=payload.wavelength_points&&payload.theta_points?payload.wavelength_points*payload.theta_points:operation.includes('kspace')?Number(payload.points||1)**2:Number(payload.points||payload.wavelength_points||1);
    const cap=Number(document.getElementById('resolutionLimit').value);
    if(!Number.isFinite(cap)||cap<9||cap>40401||count>cap)return new Response(JSON.stringify({error:`This run requests ${count} points. Adjust My maximum scan points (9–40401) or reduce the resolution.`}),{status:400});
    for(const old of [...tray.querySelectorAll('[data-finished]')].slice(0,-2))old.remove();
    const item=document.createElement('div');item.className='run-item';
    const title=document.createElement('strong');title.textContent=operation.replaceAll('_',' ');
    const cancel=document.createElement('button');cancel.textContent='Cancel';
    const bar=document.createElement('progress'),detail=document.createElement('small');
    item.append(cancel,title,bar,detail);tray.append(item);
    const started=performance.now();let jobId,finished=false;
    let estimate=null;try{estimate=Number(localStorage.getItem('timing:'+operation))*count||null;}catch{}
    const clock=setInterval(()=>{if(!jobId)detail.textContent=`Submitting · ${((performance.now()-started)/1000).toFixed(0)} s elapsed${estimate?' · prior-run estimate ~'+estimate.toFixed(0)+' s':''}`;},250);
    cancel.onclick=()=>{if(jobId){originalFetch('/api/cancel/'+jobId,{method:'POST'}).catch(()=>{});cancel.disabled=true;cancel.textContent='Stopping…';detail.textContent='Stops at the next solver checkpoint.';}};
    const abort=()=>cancel.onclick();options.signal?.addEventListener('abort',abort,{once:true});
    try{
      const response=await originalFetch(input,{...options,headers:{...options.headers,'X-Workbench-Job':'1'}});
      if(response.status!==202){finished=true;return response;}
      jobId=(await response.json()).job_id;
      for(;;){
        const poll=await originalFetch('/api/jobs/'+jobId);const job=await poll.json();
        if(!poll.ok)throw Error(job.error||'Run status unavailable');
        if(job.total){bar.max=job.total;bar.value=job.completed;}
        detail.textContent=`${job.stage} · ${job.elapsed_seconds.toFixed(0)} s elapsed${job.total?' · '+job.completed+'/'+job.total+' points':''}${job.remaining_seconds!==null?' · ~'+job.remaining_seconds.toFixed(0)+' s remaining':estimate?' · prior-run estimate ~'+estimate.toFixed(0)+' s':' · estimating after completed work'}`;
        if(['complete','failed','cancelled'].includes(job.status)){
          finished=true;item.dataset.finished='true';bar.max=1;bar.value=1;cancel.textContent='Dismiss';cancel.disabled=false;cancel.onclick=()=>item.remove();
          detail.textContent=`${job.status} · ${job.elapsed_seconds.toFixed(1)} s${job.result?.error?' · '+job.result.error:''}`;
          if(job.status==='complete')try{localStorage.setItem('timing:'+operation,String(job.elapsed_seconds/count));}catch{}
          if(job.result?.rows)setTimeout(()=>showFailedPoints(job.result.rows),100);
          if(operation==='vector_modes'&&job.status==='complete')setTimeout(()=>recordModeRun(job.result,payload),100);
          if(['tmm','spectrum'].includes(operation)&&job.status==='complete')setTimeout(()=>showTmmPhase(job.result,operation),100);
          if(operation==='dipole_ldos'&&job.status==='complete')setTimeout(()=>{const node=document.getElementById('ldosStatus');if(node)node.append(document.createTextNode(` · Collected Γ/Γ₀ = ${job.result.normalized_collected_decay_rate?.toPrecision(5)??'unavailable'}`));},100);
          return new Response(JSON.stringify(job.result),{status:job.response_status||200,headers:{'Content-Type':'application/json'}});
        }
        await new Promise(resolve=>setTimeout(resolve,650));
      }
    }catch(error){detail.textContent=error.message;cancel.textContent='Dismiss';cancel.onclick=()=>item.remove();throw error;}
    finally{clearInterval(clock);options.signal?.removeEventListener('abort',abort);if(finished&&!jobId)item.remove();refreshActions();}
  };
  function showFailedPoints(rows){
    document.getElementById('failedPointHelp')?.remove();
    const failed=rows.filter(r=>r.status&&r.status!=='converged');if(!failed.length)return;
    const panel=document.createElement('details');panel.id='failedPointHelp';panel.className='workbench-help';panel.open=true;
    const heading=document.createElement('summary');heading.textContent=`Resolve ${failed.length} failed spectral points`;panel.append(heading);
    const text=document.createElement('p');text.textContent='Open a point below, raise the Fourier budget and run Single wavelength. If a grazing cutoff is reported, shift wavelength or angle slightly. Recheck the geometry grid, then rerun the spectrum. Raw exploratory data remain available; figures require passing checks.';panel.append(text);
    for(const row of failed){const b=document.createElement('button');b.className='secondary';b.textContent=`λ ${row.wavelength_um} µm · ${row.status}`;b.onclick=()=>{document.getElementById('wavelength_um').value=row.wavelength_um;document.querySelector('.tab[data-tab="single"]').click();document.getElementById('run').scrollIntoView({block:'center'});};panel.append(b);}
    document.getElementById('quality').prepend(panel);
  }
  const legend=document.createElement('div');legend.className='material-legend';document.getElementById('structure').after(legend);
  function updateLegend(){legend.replaceChildren();const seen=new Set();for(const l of layerValues())for(const side of ['background','feature']){if(side==='feature'&&l.kind==='uniform')continue;const key=l[side+'_material'],n=l[side+'_n'],id=key==='dielectric'?key+n:key;if(seen.has(id))continue;seen.add(id);const label=document.createElement('span'),swatch=document.createElement('i');swatch.style.background=materialColor(key,n);label.append(swatch,document.createTextNode(key==='dielectric'?'Constant n='+n:key));legend.append(label);}}
  document.getElementById('layers').addEventListener('change',updateLegend);document.getElementById('preset').addEventListener('change',()=>setTimeout(updateLegend,0));updateLegend();
  const contrast=document.createElement('p');contrast.className='notice warn';contrast.hidden=true;document.getElementById('layers').after(contrast);function checkContrast(){const equal=layerValues().map((l,i)=>({l,i})).filter(({l})=>l.kind!=='uniform'&&l.background_material===l.feature_material&&(l.background_material!=='dielectric'||Number(l.background_n)===Number(l.feature_n)));contrast.hidden=!equal.length;contrast.textContent=equal.length?`Layer ${equal.map(({i})=>i+1).join(', ')} has identical background and feature material. It is optically uniform; choose different materials to create a grating.`:'';}document.getElementById('layers').addEventListener('change',checkContrast);checkContrast();
  const guideLink=document.createElement('p');guideLink.className='workbench-help';guideLink.innerHTML='<strong>Start with a worked example.</strong> <a href="/FIRST_STEPS.html" target="_blank">Open the illustrated first-use handbook</a> · literature benchmark, exact settings, expected results, and interpretation. Print or save as PDF from your browser.';
  document.querySelector('.tutorial-card')?.prepend(guideLink);
  const originalSweepPlot=drawSweep;drawSweep=function(result){originalSweepPlot(result);const c=document.getElementById('sweepPlot'),g=c.getContext('2d'),xs=result.x_values||result.rows.map(r=>r.x),ys=result.y_parameter?result.y_values:result.rows.map(r=>r.value).filter(Number.isFinite);if(!xs.length||!ys.length)return;const xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);g.fillStyle='#40576c';g.font='14px Segoe UI';for(let i=0;i<=4;i++){g.fillText((xmin+(xmax-xmin)*i/4).toPrecision(4),80+820*i/4,537);g.fillText((ymax-(ymax-ymin)*i/4).toPrecision(4),25,47+470*i/4);}c._rawData={type:'parameter_sweep',result};};
  const modeHistory=[];
  function recordModeRun(data,input){
    modeHistory.push({data,input});
    let panel=document.getElementById('modeRunHistory');if(!panel){panel=document.createElement('div');panel.id='modeRunHistory';document.getElementById('vmStatus').after(panel);}
    const table=document.createElement('table');table.innerHTML='<tr><th>Run</th><th>Actual grid</th><th>Requested mesh (µm)</th><th>Pads x/top/bottom (µm)</th><th>Mode 1 n_eff</th><th>Change</th></tr>';
    modeHistory.forEach((entry,i)=>{const row=table.insertRow();for(const text of [i+1,`${entry.data.mesh.nx} × ${entry.data.mesh.ny}`,entry.input.mesh_um,[entry.input.padding_x_um,entry.input.padding_top_um,entry.input.padding_bottom_um].join(' / '),entry.data.modes[0].n_eff.real.toFixed(7),i?(entry.data.modes[0].n_eff.real-modeHistory[i-1].data.modes[0].n_eff.real).toExponential(3):'baseline'])row.insertCell().textContent=text;});
    panel.replaceChildren(table);const note=document.createElement('p');note.className='hint';note.textContent='Compare field shapes before assigning mode identity. This history compares the first sorted eigenvalue; mesh and padding must be varied separately. Core electric fraction measures |E|² localization, not guided power.';panel.append(note);
  }
  function showTmmPhase(data,context="tmm"){
    let panel=document.getElementById(context+'PhasePanel');if(!panel){panel=document.createElement('details');panel.id=context+'PhasePanel';panel.innerHTML='<summary>Reflection and transmission phase</summary><p class="hint">Wrapped degrees of tangential electric amplitudes at the first/last interfaces; exp(−iωt). Phase is omitted at zero amplitude and for an incoherent unpolarized average.</p><canvas width="1100" height="450" class="plot"></canvas>';document.getElementById(context==='tmm'?'tmmStatus':'quality').after(panel);}
    if(data.phase_convention)panel.querySelector('p').textContent=data.phase_convention;
    const rows=data.rows||[{...data,wavelength_um:data.model.wavelength_um}],c=panel.querySelector('canvas');
    const series=['r_phase_deg','t_phase_deg'].map((key,i)=>({label:i?'Transmission phase':'Reflection phase',color:i?'#087d70':'#1768a5',x:rows.map(r=>r.wavelength_um),y:rows.map(r=>r[key]??NaN)}));
    if(!series.some(s=>s.y.some(Number.isFinite))){c.hidden=true;return;}c.hidden=false;const main=document.getElementById('plot');main.id='mainPowerPlot';c.id='plot';try{drawLines(series,'Vacuum wavelength (µm)','Wrapped phase (°)');}finally{c.id=context+'PhasePlot';main.id='plot';}c._rawData={series,phase_convention:data.phase_convention};scrollableChart(c);
  }
  // Phase requires complex amplitudes; never invent phase for R/T/A or LDOS.
  function phasePlot(canvas,values,description,axes={x0:0,x1:1,y0:0,y1:1,xlabel:'unit-cell coordinate u',ylabel:'unit-cell coordinate v',topDown:false}){
    const g=canvas.getContext('2d'),left=90,top=55,pw=canvas.width-230,ph=canvas.height-125;
    g.fillStyle='white';g.fillRect(0,0,canvas.width,canvas.height);const ny=values.length,nx=values[0].length,small=document.createElement('canvas');small.width=nx;small.height=ny;const im=small.getContext('2d').createImageData(nx,ny);
    for(let y=0;y<ny;y++)for(let x=0;x<nx;x++){const v=values[axes.topDown?y:ny-1-y][x],i=4*(y*nx+x);if(!Number.isFinite(v)){im.data.set([210,215,220,255],i);continue;}const f=Math.min(1,Math.abs(v)/180),p=Math.round(245*(1-f));im.data.set(v>=0?[230,p,p,255]:[p,p,230,255],i);}
    small.getContext('2d').putImageData(im,0,0);g.imageSmoothingEnabled=false;g.drawImage(small,left,top,pw,ph);g.strokeStyle='#234';g.strokeRect(left,top,pw,ph);g.fillStyle='#18344b';g.font='16px Segoe UI';g.fillText(description,20,28);g.fillText('−180° blue',left+pw+12,top+35);g.fillText('0° white',left+pw+12,top+65);g.fillText('+180° red',left+pw+12,top+95);g.fillText('Gray: undefined',left+pw+12,top+130);
    g.font='12px Segoe UI';for(let i=0;i<=4;i++){g.fillText((axes.x0+(axes.x1-axes.x0)*i/4).toPrecision(3),left+pw*i/4-12,top+ph+20);g.fillText((axes.y0+(axes.y1-axes.y0)*(axes.topDown?i/4:1-i/4)).toPrecision(3),left-50,top+ph*i/4+4);}g.fillText(axes.xlabel,left+pw/2-40,top+ph+42);g.save();g.translate(16,top+ph/2);g.rotate(-Math.PI/2);g.fillText(axes.ylabel,-45,0);g.restore();
    canvas._chart=null;canvas._researchHeat={values,nx,ny,left,top,pw,ph,...axes,label:description};installMapHover(canvas);canvas._rawData={quantity:description,phase_deg:values,axes};
  }
  for(const id of ['vmQuantity','fieldVariable','verticalVariable']){
    const select=document.getElementById(id);for(const name of ['Ex','Ey','Ez','Hx','Hy','Hz'])select.add(new Option(name+' phase (°)',name+'_phase_deg'));
    const previous=select.onchange;select.onchange=function(){if(!this.value.endsWith('_phase_deg'))return previous?.call(this);const values=id==='vmQuantity'?vectorModeData?.modes[Number(document.getElementById('vmMode').value)||0]?.fields[this.value]:id==='fieldVariable'?fieldData?.[this.value]:verticalData?.[this.value];if(values){let axes;if(id==='vmQuantity')axes={x0:vectorModeData.x_um[0],x1:vectorModeData.x_um.at(-1),y0:vectorModeData.y_um[0],y1:vectorModeData.y_um.at(-1),xlabel:'x (µm)',ylabel:'y (µm)'};if(id==='verticalVariable')axes={x0:0,x1:1,y0:verticalData.z_um[0],y1:verticalData.z_um.at(-1),xlabel:'position along primitive vector',ylabel:'z (µm)',topDown:true};const canvas=document.getElementById(id==='vmQuantity'?'vmPlot':'plot');phasePlot(canvas,values,this.selectedOptions[0].text,axes);if(id==='verticalVariable'){const h=canvas._researchHeat;drawVerticalStructureProfile(canvas.getContext('2d'),verticalData,h.left,h.top,h.pw,h.ph);}}};
  }
  const polSelect=document.getElementById('polMapQuantity');for(const name of ['Ep','Es'])polSelect.add(new Option(name+' phase (°)',name+'_phase_deg'));const oldPol=polSelect.onchange;polSelect.onchange=function(){if(!this.value.endsWith('_phase_deg'))return oldPol?.call(this);if(lastPolarizationMap){const d=lastPolarizationMap;phasePlot(document.getElementById('polarizationMapCanvas'),d.maps[this.value],this.selectedOptions[0].text+' · outgoing local p/s basis',{x0:d.u[0],x1:d.u.at(-1),y0:d.v[0],y1:d.v.at(-1),xlabel:'kx / (nᵢ k₀)',ylabel:'ky / (nᵢ k₀)'});}};
  const mapSelect=document.getElementById('mapQuantity');mapSelect.add(new Option('Specular reflection phase (°)','r_phase_deg'));mapSelect.add(new Option('Specular transmission phase (°)','t_phase_deg'));const oldScalar=drawScalarMap;drawScalarMap=function(data){if(!data.quantity.endsWith('_phase_deg'))return oldScalar(data);lastScatteringMap=data;const angle=data.kind==='angle_wavelength',x=angle?data.wavelength_um:data.u,y=angle?data.theta_deg:data.v;phasePlot(document.getElementById('scatteringMapCanvas'),data.values,data.quantity+' · coherent component; check phase convergence separately',{x0:x[0],x1:x.at(-1),y0:y[0],y1:y.at(-1),xlabel:angle?'Vacuum wavelength (µm)':'kx / (nᵢ k₀)',ylabel:angle?'Incidence angle (°)':'ky / (nᵢ k₀)'});};
  const branch=document.getElementById('runQBranches');if(branch){const link=document.createElement('button');link.type='button';link.className='secondary';link.textContent='Track several resonance branches →';link.onclick=()=>{showWorkspace('sweeps');branch.closest('details').open=true;branch.scrollIntoView({block:'center'});refreshActions();};document.getElementById('runMultiResonance')?.parentElement.append(link);if(!link.isConnected)document.getElementById('runResonance').parentElement.append(link);}
  const vectorRun=document.getElementById('runVectorMode');const vectorShape=document.createElement('div');vectorShape.className='row';vectorShape.innerHTML='<div><label for="vmCoreShape">Core cross-section</label><select id="vmCoreShape"><option value="rectangle">Rectangle</option><option value="trapezoid">Trapezoid</option><option value="ellipse">Ellipse</option></select></div><div><label for="vmSidewall">Sidewall angle (°)</label><input id="vmSidewall" type="number" min="20" max="90" value="90"></div>';vectorRun.closest('.card').querySelector('.actions')?.before(vectorShape);
  const refine=document.createElement('button');refine.id='runModeRefinement';refine.className='secondary';refine.textContent='Run mesh + padding comparison';vectorRun.after(refine);
  refine.onclick=async()=>{refine.disabled=true;const read=id=>Number(document.getElementById(id).value),input={wavelength_um:read('vmLambda'),core_width_um:read('vmWidth'),core_height_um:read('vmHeight'),core_n:read('vmCore'),substrate_n:read('vmSub'),cladding_n:read('vmClad'),padding_x_um:read('vmPadX'),padding_top_um:read('vmPadTop'),padding_bottom_um:read('vmPadBottom'),mesh_um:read('vmMesh'),modes:read('vmModes'),guess:read('vmGuess'),boundary:document.getElementById('vmBoundary').value,core_shape:document.getElementById('vmCoreShape').value,sidewall_angle_deg:read('vmSidewall')};try{for(const test of [input,{...input,mesh_um:input.mesh_um/2},{...input,padding_x_um:input.padding_x_um*1.25,padding_top_um:input.padding_top_um*1.25,padding_bottom_um:input.padding_bottom_um*1.25}]){const response=await fetch('/api/vector_modes',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(test)});const result=await response.json();if(!response.ok)throw Error(result.error);vectorModeData=result;document.getElementById('vmMode').replaceChildren(...result.modes.map((m,i)=>new Option('Mode '+(i+1),i)));drawVectorMode();}document.getElementById('vmStatus').textContent='Three runs completed: baseline, half mesh at fixed padding, and larger padding at baseline mesh. Inspect mode identity and numerical changes below.';}catch(e){document.getElementById('vmStatus').textContent=e.message;}finally{refine.disabled=false;}};

  // Resonance evidence uses the same complex field payload as the vertical-field workspace.
  const resCase=document.getElementById('resFieldCase');
  if(resCase){
    const quantity=document.createElement('select');quantity.id='resFieldQuantity';
    for(const [value,label] of [['E2','Electric intensity |E|²'],['H2','Magnetic intensity |H|²'],['Ex_phase_deg','Eₓ phase (°)'],['Ey_phase_deg','Eᵧ phase (°)'],['Ez_phase_deg','Eᶻ phase (°)'],['Hx_phase_deg','Hₓ phase (°)'],['Hy_phase_deg','Hᵧ phase (°)'],['Hz_phase_deg','Hᶻ phase (°)']])quantity.add(new Option(label,value));
    const wrap=document.createElement('div');wrap.innerHTML='<label for="resFieldQuantity">Displayed resonance field</label>';wrap.append(quantity);resCase.parentElement.after(wrap);
    quantity.onchange=()=>{if(!resonanceFieldData)return;const item=resonanceFieldData.cases[Number(resCase.value)||0],key=quantity.value;if(key==='E2')return drawResonanceFieldCase();const values=item.field[key];if(!values)return;phasePlot(document.getElementById('resFieldCanvas'),values,quantity.selectedOptions[0].text+' · '+item.label,{x0:0,x1:1,y0:item.field.z_um[0],y1:item.field.z_um.at(-1),xlabel:item.field.plane==='xz'?'a₁ coordinate':'a₂ coordinate',ylabel:'z (µm; top interface = 0)',topDown:true});const h=document.getElementById('resFieldCanvas')._researchHeat;drawVerticalStructureProfile(document.getElementById('resFieldCanvas').getContext('2d'),item.field,h.left,h.top,h.pw,h.ph);};
    resCase.addEventListener('change',()=>quantity.dispatchEvent(new Event('change')));
  }

  // Compact long workspaces without hiding the calculation title.
  const compact=document.createElement('button');compact.type='button';compact.className='secondary';compact.textContent='Collapse result cards';
  compact.onclick=()=>{const cards=[...document.querySelectorAll('section .card')].filter(card=>!card.hidden);const collapse=cards.some(card=>!card.classList.contains('compact-collapsed'));for(const card of cards)card.classList.toggle('compact-collapsed',collapse);compact.textContent=collapse?'Expand result cards':'Collapse result cards';};
  nav.append(compact);
  for(const card of document.querySelectorAll('section .card')){const h=card.querySelector(':scope > h2');if(!h)continue;const toggle=document.createElement('button');toggle.type='button';toggle.className='card-toggle';toggle.title='Collapse or expand this card';toggle.textContent='▾';toggle.onclick=()=>{card.classList.toggle('compact-collapsed');toggle.textContent=card.classList.contains('compact-collapsed')?'▸':'▾';};h.prepend(toggle);}

  // Searchable control glossary assembled from visible labels plus curated physics guidance.
  const definitions={
    'Requested Fourier harmonics':'Controls the lateral Fourier basis. Increase until the observable and field pattern stop changing. Runtime and memory grow rapidly.',
    'Geometry grid per axis':'Samples material boundaries before Fourier factorization. Refine independently from the harmonic budget.',
    'Incident angle θ (°)':'Angle from the incident-side surface normal. Avoid an exact diffraction-order cutoff when convergence warnings report grazing propagation.',
    'Azimuth φ (°)':'Direction of the in-plane incident wavevector relative to primitive vector a₁.',
    'Points':'Number of actually solved samples. Display smoothing does not add physical resolution.',
    'Mesh step (µm)':'Finite-difference cell size. Repeat at a smaller value while holding the domain padding fixed.',
    'Maximum scan points':'Personal guardrail for total solved points. Raise it only after a coarse scan estimates cost and locates the feature.',
    'Noise σ':'Standard deviation assumed by the independent Gaussian likelihood. Estimate it from measurement repeats or a justified noise model.',
    'Bootstrap replicates':'Repeated residual resamples used to expose fit instability. Use at least 100 for a final interval when runtime allows.',
    'Refinement rounds':'Number of interval-narrowing passes. More rounds help only when every round still samples the linewidth.',
    'Points per round':'Spectral samples in each adaptive pass. A fitted linewidth needs several calculated points across its width.',
    'Odd grid points':'Centered k-space sampling. Odd counts include the origin exactly.',
    'Minimum channel power':'Masks normalized polarization where division by nearly zero outgoing power is unstable.',
    'Vertical grid points':'Samples the one-dimensional mode eigenproblem. Refine until n_eff and the field profile agree.',
    'Maximum modes':'Number of eigenpairs requested near the supplied index guess; extra modes increase solve cost and help identify crossings.',
    'NA':'Objective numerical aperture in the upper medium. The accepted half-angle is asin(NA/n).'
  };
  const help=document.createElement('dialog');help.id='parameterHelp';help.innerHTML='<form method="dialog" class="help-head"><h2>Search controls and numerical guidance</h2><button aria-label="Close">Close</button></form><input id="helpSearch" type="search" placeholder="Search wavelength, Fourier, mesh, Q, NA…"><div id="helpResults"></div><p><a href="/FIRST_STEPS.html" target="_blank">Open the complete worked guide ↗</a></p>';document.body.append(help);
  const helpButton=document.createElement('button');helpButton.type='button';helpButton.className='secondary';helpButton.textContent='Search help';nav.append(helpButton);
  function helpEntries(){const seen=new Set(),items=[];for(const label of document.querySelectorAll('label')){const name=label.textContent.trim();if(!name||seen.has(name))continue;seen.add(name);const control=label.htmlFor?document.getElementById(label.htmlFor):label.parentElement?.querySelector('input,select,textarea');const card=label.closest('.card');const range=control?.type==='number'?`Range ${control.min||'unbounded'} to ${control.max||'unbounded'}${control.step?' · step '+control.step:''}.`:'';const local=label.parentElement?.querySelector('.hint')?.textContent||'';items.push({name,section:card?.querySelector('h2')?.textContent.trim()||'General',text:definitions[name]||local||'Use the displayed units and repeat the calculation at higher numerical resolution before quoting a sensitive result.',range});}return items;}
  function renderHelp(){const q=document.getElementById('helpSearch').value.trim().toLowerCase(),items=helpEntries().filter(v=>(v.name+' '+v.section+' '+v.text).toLowerCase().includes(q)).slice(0,80),out=document.getElementById('helpResults');out.innerHTML=items.length?items.map(v=>`<article><h3>${v.name}</h3><small>${v.section}</small><p>${v.text} ${v.range}</p></article>`).join(''):'<p>No matching control. Try a physical quantity or solver name.</p>';}
  helpButton.onclick=()=>{renderHelp();help.showModal();setTimeout(()=>document.getElementById('helpSearch').focus(),0);};document.getElementById('helpSearch').oninput=renderHelp;

  // Durable run history survives a server restart; numerical result arrays remain in project exports.
  const historyButton=document.createElement('button');historyButton.type='button';historyButton.className='secondary';historyButton.textContent='Run history';nav.append(historyButton);
  const historyDialog=document.createElement('dialog');historyDialog.innerHTML='<form method="dialog" class="help-head"><h2>Recent calculations</h2><button>Close</button></form><div id="jobHistory">Loading…</div>';document.body.append(historyDialog);
  historyButton.onclick=async()=>{historyDialog.showModal();const out=document.getElementById('jobHistory');try{const response=await originalFetch('/api/job-history'),data=await response.json();out.innerHTML=data.jobs.length?'<table><tr><th>Operation</th><th>Status</th><th>Finished</th><th>Elapsed</th><th>Progress</th></tr>'+data.jobs.map(j=>`<tr><td>${j.operation||'calculation'}</td><td>${j.status}</td><td>${j.finished_at?new Date(j.finished_at*1000).toLocaleString():'active'}</td><td>${Number(j.elapsed_seconds||0).toFixed(1)} s</td><td>${j.total?j.completed+'/'+j.total:j.stage||'—'}</td></tr>`).join('')+'</table>':'<p>No recorded calculations yet.</p>';}catch(error){out.textContent=error.message;}};
  const bayesSigma=document.getElementById('bayesSigma');if(bayesSigma){const correlation=document.createElement('div');correlation.innerHTML='<label for="bayesCorrelation">AR(1) wavelength-noise correlation</label><input id="bayesCorrelation" type="number" min="-0.94" max="0.94" step="0.01" value="0"><p class="hint">Use zero for independent errors. Estimate nonzero correlation from repeated-measurement residuals rather than tuning it for a narrower interval.</p>';bayesSigma.parentElement.after(correlation);}
  const optimizerSeed=document.getElementById('optimizerSeed');if(optimizerSeed){const robust=document.createElement('div');robust.className='row';robust.innerHTML='<div><label for="optimizerRobustSamples">Fabrication samples / candidate</label><input id="optimizerRobustSamples" type="number" min="1" max="12" value="1"></div><div><label for="optimizerRobustWeight">Variability penalty</label><input id="optimizerRobustWeight" type="number" min="0" max="2" step="0.05" value="0.25"></div>';optimizerSeed.closest('.row3').after(robust);}
  const bandCard=document.getElementById('bandModeCanvas')?.closest('.card');if(bandCard){const zone=document.createElement('details');zone.innerHTML='<summary>Complete-gap screen over the full reciprocal cell</summary><p class="hint">The standard Γ–X–M–Γ path can miss a band extremum. This calculation samples one complete reciprocal primitive cell before reporting a candidate complete gap.</p><div class="row"><div><label for="fullZonePoints">Odd k points per reciprocal axis</label><input id="fullZonePoints" type="number" min="5" max="31" step="2" value="9"></div><div class="actions"><button id="runFullZone" class="primary">Screen complete gaps</button></div></div><p id="fullZoneStatus" class="notice">Start with 9 × 9, then repeat at 15 × 15 or finer.</p><div id="fullZoneReport"></div>';bandCard.append(zone);document.getElementById('runFullZone').onclick=async()=>{const button=document.getElementById('runFullZone'),note=document.getElementById('fullZoneStatus');button.disabled=true;note.textContent='Scanning the reciprocal primitive cell…';try{const response=await fetch('/api/full_zone_gaps',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({model:bandModel(),points:Number(document.getElementById('fullZonePoints').value)})}),data=await response.json();if(!response.ok)throw Error(data.error);note.className='notice '+(data.complete_sampled_gaps.length?'ok':'warn');note.textContent=`Solved ${data.grid_points_per_axis} × ${data.grid_points_per_axis} k points. ${data.complete_sampled_gaps.length} sampled complete polarization gap(s).`;document.getElementById('fullZoneReport').innerHTML=data.complete_sampled_gaps.length?'<table><tr><th>Polarization</th><th>Bands</th><th>a/λ interval</th><th>Relative width</th><th>λ interval (µm)</th></tr>'+data.complete_sampled_gaps.map(g=>`<tr><td>${g.polarization}</td><td>${g.between_bands.join('–')}</td><td>${g.lower_a_over_lambda.toFixed(5)}–${g.upper_a_over_lambda.toFixed(5)}</td><td>${(100*g.relative_width).toFixed(2)}%</td><td>${g.wavelength_um.map(v=>v.toFixed(5)).join('–')}</td></tr>`).join('')+'</table><p>'+data.warning+'</p>':'<p>No complete gap was found on this grid. A directional path gap may still exist.</p>';}catch(error){note.className='notice bad';note.textContent=error.message;}finally{button.disabled=false;}};}
  refreshActions();
})();
