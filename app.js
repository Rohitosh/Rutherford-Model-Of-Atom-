  // app.js
  // Place in same folder as index.html and styles.css

  (() => {
    // constants
    const e_charge = 1.602176634e-19;
    const epsilon0 = 8.8541878128e-12;
    const pi = Math.PI;
    const k_coulomb = 1.0 / (4.0 * pi * epsilon0);

    // DOM
    const simCanvas = document.getElementById('simCanvas');
    const ctx = simCanvas.getContext('2d');
    const histCanvas = document.getElementById('histCanvas');
    const hctx = histCanvas.getContext('2d');
    const btnSim = document.getElementById('btnSim');
    const btnPlay = document.getElementById('btnPlay');
    const btnPause = document.getElementById('btnPause');
    const btnReset = document.getElementById('btnReset');
    const btnExport = document.getElementById('btnExport');
    const inputN = document.getElementById('inputN');
    const inputE = document.getElementById('inputE');
    const inputB = document.getElementById('inputB');
    const inputBins = document.getElementById('inputBins');
    const statParticles = document.getElementById('statParticles');
    const speedRange = document.getElementById('speed');

    let particles = []; // particle objects for animation
    let simState = { playing: false, frame:0, speed:1.0 };
    let simParams = {};

    // canvas geometry mapping: we will place foil plane at center X
    const W = simCanvas.width;
    const H = simCanvas.height;
    const foilXpx = Math.round(W * 0.5);
    const foilLen = 220;
    const nucleusRadiusPx = 8;
    const detectorRpx = 220;

    // utility: map physical meters to pixels - we choose a scale based on bmax and display region
    function computeScale(bmax_m) {
      // we want bmax to fit within ~ half canvas height
      const pixelsPerMeter = (H * 0.9) / (2.0 * bmax_m); // rough
      return pixelsPerMeter;
    }

    // run the analytic Monte Carlo simulation and produce particle trajectories and histogram
    function runSimulation() {
      const N = Math.max(100, parseInt(inputN.value, 10));
      const E_MeV = parseFloat(inputE.value);
      const bmax_ang = parseFloat(inputB.value); // Å units in UI
      const bins = parseInt(inputBins.value, 10);

      const bmax = bmax_ang * 1e-10; // convert Å -> meters
      const E_joule = E_MeV * 1e6 * e_charge;
      // choose arbitrary non-relativistic v for animation (scales only)
      const v0 = Math.sqrt(2.0 * E_joule / (4.0 * 1.66053906660e-27)); // alpha mass 4 amu

      // theoretical prefactor: (k q1 q2 / (16 E))^2 but we will compute later for dσ/dΩ normalized
      const Z1 = 2.0, Z2 = 79.0;
      const q1 = Z1 * e_charge, q2 = Z2 * e_charge;
      const pref = Math.pow((k_coulomb * q1 * q2) / (4.0 * E_joule), 2.0); // note: using 1/(4E) to match formula; will scale later

      // prepare histogram
      const hist = new Array(bins).fill(0);
      const thetaCenters = new Array(bins);
      for (let i=0;i<bins;i++) thetaCenters[i] = (i + 0.5) * Math.PI / bins; // radians

      // scale: meters->pixels
      const scale = computeScale(bmax);
      const pxPerM = scale;

      // produce N particles
      particles = [];
      const angles_deg = [];
      for (let i=0;i<N;i++){
        // sample b with area-weighting
        const u = Math.random();
        const b = bmax * Math.sqrt(u);
        // sign
        const sign = (Math.random() < 0.5) ? 1.0 : -1.0;
        const y0 = sign * b;

        // compute scattering angle analytically
        const t2 = (k_coulomb * q1 * q2) / (2.0 * E_joule * b);
        let theta = 2.0 * Math.atan(t2);
        if (!isFinite(theta)) theta = Math.PI;
        angles_deg.push(theta * 180.0 / Math.PI);

        // histogram bin
        let deg = theta * 180.0 / Math.PI;
        let bin = Math.floor(deg / (180.0 / bins));
        if (bin < 0) bin = 0;
        if (bin >= bins) bin = bins-1;
        hist[bin]++;

        // create a trajectory: we will record a simple pre-foil line then post-foil line at angle theta
        // position units are in meters; we'll map to pixels in drawing
        // start x at -Xstart (m)
        const startX = - (W/2) / pxPerM * 0.8; // map half-canvas to meters left
        const foilX = 0.0;
        // pre-positions: many frames approaching foil
        const preFrames = 24;
        const postFrames = 120;
        const traj = [];
        const dxPre = (foilX - startX) / preFrames;
        for (let f=0; f<preFrames; f++){
          const x = startX + dxPre * f;
          traj.push({x:x, y:y0});
        }
        // at foil: compute post velocity direction
        const thetaSigned = (y0 >= 0 ? theta : -theta);
        const vx = v0 * Math.cos(thetaSigned);
        const vy = v0 * Math.sin(thetaSigned);
        // step in time for post frames
        let x = foilX, y = y0;
        const dt = 1e-17; // abstract time per frame
        for (let f=0; f<postFrames; f++){
          x += vx * dt;
          y += vy * dt;
          traj.push({x:x, y:y});
          // allow break if left region
          if (Math.abs(x) > 1e-12) {
            // allow continue — still keep trajectory
          }
        }
        particles.push({traj:traj, color: (Math.abs(thetaSigned) > 0.2 ? '#ff6b6b' : '#6eb5ff')});
      }

      // update UI
      statParticles.textContent = N;

      // compute theoretical pdf per theta bin: dσ/dΩ ~ 1/sin^4(θ/2). We'll convert to expected counts per bin
      const theory = new Array(bins).fill(0.0);
      let norm = 0.0;
      for (let i=0;i<bins;i++){
        const theta = thetaCenters[i];
        const sinHalf = Math.sin(theta/2.0);
        if (sinHalf <= 0) { theory[i]=0; continue;}
        const dsdo = pref / Math.pow(sinHalf, 4.0); // unnormalized
        // convert to counts per theta by multiplying with 2π sinθ dθ
        const dtheta = Math.PI / bins;
        const mass = dsdo * (2.0 * Math.PI * Math.sin(theta)) * dtheta;
        theory[i] = mass;
        norm += mass;
      }
      const expected = theory.map(v => (norm>0 ? v / norm * N : 0));

      // save sim results to state
      simParams = {N, E_MeV, bmax, bins, hist, expected, thetaCenters, angles:angles_deg, particles, pxPerM};

      // draw histogram and prepare animation
      drawHistogram(simParams);
      simState.frame = 0;
    }

    // draw histogram: histogram (black background), simulated bars vs theoretical curve
    function drawHistogram(state) {
      const w = histCanvas.width, h = histCanvas.height;
      hctx.clearRect(0,0,w,h);
      // background
      hctx.fillStyle = '#070707';
      hctx.fillRect(0,0,w,h);

      const bins = state.bins;
      const maxCount = Math.max(...state.hist, 1);
      const barW = w / bins;
      // draw simulated bars
      for (let i=0;i<bins;i++){
        const barH = (state.hist[i] / maxCount) * (h - 30);
        hctx.fillStyle = '#000000';
        hctx.fillRect(i*barW, h-20-barH, barW*0.9, barH);
      }
      // draw theoretical line scaled to same max
      const maxExp = Math.max(...state.expected,1);
      hctx.beginPath();
      hctx.strokeStyle = '#ff8b3a';
      hctx.lineWidth = 2;
      for (let i=0;i<bins;i++){
        const x = i*barW + barW*0.5;
        const y = h - 20 - (state.expected[i]/maxExp)*(h-30);
        if (i==0) hctx.moveTo(x,y); else hctx.lineTo(x,y);
      }
      hctx.stroke();

      // axes
      hctx.fillStyle = '#ddd';
      hctx.font = '12px monospace';
      hctx.fillText('θ (deg)', w/2 - 20, h-2);
    }

    // animate the experiment canvas
    function animateSim() {
      ctx.clearRect(0,0,W,H);
      // background
      ctx.fillStyle = '#060606';
      ctx.fillRect(0,0,W,H);

      // draw foil (vertical line at center)
      ctx.strokeStyle = 'gold';
      ctx.lineWidth = 4;
      ctx.beginPath();
      ctx.moveTo(foilXpx, H/2 - foilLen/2);
      ctx.lineTo(foilXpx, H/2 + foilLen/2);
      ctx.stroke();

      // draw nucleus at center
      ctx.beginPath();
      ctx.fillStyle = '#ffcf4d';
      ctx.arc(foilXpx, H/2, nucleusRadiusPx, 0, 2*pi);
      ctx.fill();
      ctx.strokeStyle = '#ffaa00';
      ctx.lineWidth = 1;
      ctx.stroke();

      // detector circle
      ctx.strokeStyle = 'rgba(110,181,255,0.25)';
      ctx.lineWidth = 2;
      ctx.beginPath();
      ctx.arc(foilXpx, H/2, detectorRpx, 0, 2*pi);
      ctx.stroke();

      // source text
      ctx.fillStyle = '#ddd';
      ctx.font = '14px sans-serif';
      ctx.fillText('Alpha source (left)', 20, 26);

      // draw each particle at its current frame
      const speed = parseFloat(speedRange.value);
      simState.speed = speed;

      if (simParams.particles && Array.isArray(simParams.particles)) {
        for (let pi=0; pi<simParams.particles.length; pi++){
          const traj = simParams.particles[pi].traj;
          const len = traj.length;
          // frame index scaled by speed
          let idx = Math.floor(simState.frame * (1 + speed*0.2));
          if (idx >= len) idx = len-1;
          const pos = traj[idx];
          // map meters to screen
          const sx = foilXpx + pos.x * simParams.pxPerM;
          const sy = H/2 - pos.y * simParams.pxPerM;
          // draw trail (few previous points)
          ctx.beginPath();
          ctx.strokeStyle = simParams.particles[pi].color;
          ctx.lineWidth = 1.6;
          for (let t = Math.max(0, idx-6); t<=idx; t++){
            const p = traj[t];
            const tx = foilXpx + p.x * simParams.pxPerM;
            const ty = H/2 - p.y * simParams.pxPerM;
            if (t===Math.max(0, idx-6)) ctx.moveTo(tx,ty); else ctx.lineTo(tx,ty);
          }
          ctx.stroke();

          // dot
          ctx.fillStyle = simParams.particles[pi].color;
          ctx.beginPath();
          ctx.arc(sx, sy, 3, 0, 2*pi);
          ctx.fill();
        }
      }

      // UI overlay
      ctx.fillStyle = '#ddd';
      ctx.font = '12px monospace';
      ctx.fillText(`frame: ${simState.frame}`, 12, H-8);

      if (simState.playing) {
        simState.frame++;
        // wrap
        if (simState.frame > 6000) simState.frame = 0;
      }
      requestAnimationFrame(animateSim);
    }

    // export CSV of angles (for external plotting)
    function exportAnglesCSV() {
      if (!simParams.angles) return;
      let csv = 'particle,theta_deg\n';
      simParams.angles.forEach((v,i) => csv += `${i},${v}\n`);
      const url = URL.createObjectURL(new Blob([csv], {type: 'text/csv'}));
      const a = document.createElement('a');
      a.href = url; a.download = 'angles_export.csv';
      document.body.appendChild(a); a.click(); a.remove();
    }

    // wire up UI
    btnSim.addEventListener('click', () => { runSimulation(); });
    btnPlay.addEventListener('click', () => { simState.playing = true; });
    btnPause.addEventListener('click', () => { simState.playing = false; });
    btnReset.addEventListener('click', () => { simState.playing = false; simState.frame = 0; });
    btnExport.addEventListener('click', exportAnglesCSV);

    // start animation loop
    requestAnimationFrame(animateSim);
  })();
