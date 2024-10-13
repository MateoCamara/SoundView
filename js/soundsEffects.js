const audioContext = new (window.AudioContext || window.webkitAudioContext)();
let gainNode;
let convolverNode;
let audioBuffer;
let audioSource;

speedControl.addEventListener('input', function() {
    const speed = parseFloat(this.value).toFixed(2);
    speedValue.textContent = `x${speed}`;
    if (audioSource) {
        audioSource.playbackRate.value = this.value;
    }
});

gainControl.addEventListener('input', function() {
    const volume = parseFloat(this.value);
    gainValue.textContent = volume.toFixed(1);
    if (gainNode) {
        gainNode.gain.value = volume;
    }
});

function playAudio() {
    audio_manager.stopAllBufferNodes();
    if (audioSource) {
        audioSource.stop();
    }
    
    fetch(selectedSound.preview_url)
    .then(response => response.arrayBuffer())
    .then(arrayBuffer => audioContext.decodeAudioData(arrayBuffer))
    .then(decodedAudio => {
        audioBuffer = decodedAudio;

        audioSource = audioContext.createBufferSource();
        audioSource.buffer = audioBuffer;
        audioSource.playbackRate.value = speedControl.value;

        // Crear y conectar el nodo de ganancia
        gainNode = audioContext.createGain();
        gainNode.gain.value = gainControl.value;
        audioSource.connect(gainNode);
        gainNode.connect(audioContext.destination);

        audioSource.onended = function() {
            updateProgress(1);
        };

        audioSource.start();

        // Progreso de reproducción en tiempo real
        let animationId = requestAnimationFrame(updateProgressLoop);

        function updateProgressLoop() {
            const progress = audioContext.currentTime / audioBuffer.duration;
            updateProgress(progress);

            if (progress >= 1) {
                cancelAnimationFrame(animationId);
            } else {
                animationId = requestAnimationFrame(updateProgressLoop);
            }
        }
    })
    .catch(error => console.error('Error loading audio:', error));
}
