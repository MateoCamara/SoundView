function applyEffects(index, value) {
  const audioContext = new AudioContext();
  const sourceNode = audioContext.createBufferSource();
  const gainNode = audioContext.createGain();
  const filterNode = audioContext.createBiquadFilter();
  const analyserNode = audioContext.createAnalyser();

  sourceNode.connect(gainNode);
  gainNode.connect(filterNode);
  filterNode.connect(analyserNode);
  analyserNode.connect(audioContext.destination);

  const audioBuffer = audioContext.createBuffer(1, currentSound.length, audioContext.sampleRate);
  const channelData = audioBuffer.getChannelData(0);
  channelData.set(currentSound);

  sourceNode.buffer = audioBuffer;
  sourceNode.start();

  const initialData = new Float32Array(analyserNode.fftSize);
  analyserNode.getFloatTimeDomainData(initialData);

  const areInputsAtMidValue = Array.from(transformInputs).every(
    (input) => parseFloat(input.value) === 0
  );

  if (areInputsAtMidValue) {
    sourceNode.playbackRate.value = 1;
    filterNode.type = 'allpass';
    gainNode.disconnect();
    gainNode.connect(filterNode);
  } else {
      switch (index) {
        case 0:
          // Aplicar efecto de velocidad
          sourceNode.playbackRate.value = value;
          break;
        case 1:
          // Aplicar efecto de ecualización
          filterNode.type = 'peaking';
          filterNode.frequency.value = value;
          break;
        case 2:
          // Aplicar efecto de distorsión
          if (value > 0) {
            const distortionNode = audioContext.createWaveShaper();
            gainNode.disconnect();
            gainNode.connect(distortionNode);
            distortionNode.connect(filterNode);
          } else {
            gainNode.disconnect();
            gainNode.connect(filterNode);
          }
          break;
        case 3:
          // Aplicar efecto de inversión
          if (value > 0) {
            const invertNode = audioContext.createGain();
            invertNode.gain.value = -1;
            gainNode.disconnect();
            gainNode.connect(invertNode);
            invertNode.connect(filterNode);
          } else {
            gainNode.disconnect();
            gainNode.connect(filterNode);
          }
          break;
      }
    }
}