function adjustAudioToExpectedSize(
  waveform,
  sampleRate = 22050 * 2,
  duration = 0.742
) {
  let expectedDuration = math.floor(sampleRate * duration);
  let adjustedWaveform = [];
  if (waveform.length >= expectedDuration) {
    adjustedWaveform = waveform.slice(0, expectedDuration);
  } else {
    let auxZeros = new Array(expectedDuration - waveform.length).fill(0);
    adjustedWaveform = waveform.push(...auxZeros);
  }
  return adjustedWaveform;
}

function getWaveformFromPreview(callback, audioPreviewUrl, sampleRate = 22050) {
  let wavesurfer = WaveSurfer.create({
    container: "#waveform",
    scrollParent: true,
  });

  wavesurfer.load(audioPreviewUrl);
  wavesurfer.on(
    "ready",
    () =>
      wavesurfer
        .exportPCM(sampleRate, 10000, true)
        .then((waveform) => callback(waveform)),
    sampleRate
  );
}

function encodeAudio(arangedSpectrogram) {
  let batchInputShape = encoderModel.inputLayers[0].batchInputShape; // In the default model, this is 512 x 64
  let tensor = tf.tensor2d(
    arangedSpectrogram,
    [batchInputShape[1], batchInputShape[2]],
    "float32"
  );
  let [mu_tensor, log_variance_tensor] = encoderModel.predict(
    tf.reshape(tensor, (shape = [1, batchInputShape[1], batchInputShape[2], 1]))
  );
  let mu_latent_space = mu_tensor.dataSync();
  let log_variance_latent_space = log_variance_tensor.dataSync();
  return { mu_latent_space, log_variance_latent_space };
}

async function getUserUploadedModel() {
  uploadVAEs.style.display = "none";
  const uploadEncoderJson = document.getElementById("upload-encoder-json");
  const uploadEncoderWeights = document.getElementById(
    "upload-encoder-weights"
  );
  const uploadDecoderJson = document.getElementById("upload-decoder-json");
  const uploadDecoderWeights = document.getElementById(
    "upload-decoder-weights"
  );
  if (
    !uploadEncoderJson ||
    !uploadEncoderWeights ||
    !uploadDecoderJson ||
    !uploadDecoderWeights
  ) {
    logInfo("Failed to include user VAE!!");
    return;
  }
  console.log(uploadDecoderWeights);
  encoderModel = await tf.loadLayersModel(
    tf.io.browserFiles([
      uploadEncoderJson.files[0],
      ...uploadEncoderWeights.files,
    ])
  );
  logInfo("Loaded user encoder.");
  decoderModel = await tf.loadLayersModel(
    tf.io.browserFiles([
      uploadDecoderJson.files[0],
      ...uploadDecoderWeights.files,
    ])
  );
  logInfo("Loaded user decoder.");
  start();
}

function phmod(ph) {
  return ph < 0 ? PI2 + (ph % PI2) : ph % PI2;
}

function unwrap(input, output) {
  const size = input.length;
  if (!output) output = input;
  if (output === true) output = new Array(size);

  let shift = 0;
  let prev = phmod(input[0]);
  output[0] = prev;
  for (let i = 1; i < size; i++) {
    const current = phmod(input[i]);
    const diff = current - prev;
    if (diff < -PI) shift += PI2;
    else if (diff > PI) shift -= PI2;
    output[i] = current + shift;
    prev = current;
  }
  return output;
}

function normalizeSpec(data, minVal, maxVal, minData, maxData) {
  let normArray = math.dotDivide(
    math.subtract(data, minData),
    maxData - minData
  );
  normArray = math.add(math.dotMultiply(normArray, maxVal - minVal), minVal);
  return normArray;
}

function denormalizeSpec(normData, minVal, maxVal, minData, maxData) {
  let array = math.dotDivide(
    math.subtract(normData, minVal),
    math.subtract(maxVal, minVal)
  );
  array = math.add(
    math.dotMultiply(array, math.subtract(maxData, minData)),
    minData
  );
  return array;
}

// function magphase(complexValue) {
//     let angle = complexValue.arg();
//let unwrappedPhase = unwrap(angle, true);
//     let magnitude = complexValue.abs();
//     magnitude <= 0 ? magnitude = 0.00001: magnitude // check if it's zero or negative (convert to -100 db)
//     let dbMag = 20 * math.log10(magnitude);
//     return [dbMag, angle];
// }

function spectrogram(spectrogram, type = "2D") {
  let spectrogramResult = [];
  let magSpec = [];
  if (type === "2D") {
    for (let frame = 0; frame < spectrogram.length; frame++) {
      let magSpecAux = [];
      for (let index = 0; index < spectrogram[frame].length; index++) {
        let magnitude = spectrogram[frame][index];
        magnitude <= 0 ? (magnitude = 0.00001) : magnitude; // check if it's zero or negative (convert to -100 db)
        let mag = 20 * math.log10(magnitude);
        magSpecAux.push(mag);
      }
      magSpec.push(magSpecAux);
    }
    let magSpecNorm = normalizeSpec(magSpec, 0, 1, -100, 0);
    spectrogramResult = magSpecNorm;
  }
  return spectrogramResult;
}

async function convertPredictedSpectrogramIntoAudio(
  callback,
  predictedSpec,
  type = "2D",
  framelength = 1024
) {
  let magSpecNorm = predictedSpec.slice(0, predictedSpec.length / 2);
  let phaSpecNorm = predictedSpec.slice(
    predictedSpec.length / 2,
    predictedSpec.length
  );
  let magSpec = denormalizeSpec(magSpecNorm, 0, 1, -100, 0);
  let phaSpec = denormalizeSpec(phaSpecNorm, 0, 1, -100, 100);

  let reconstructedSpec = [];
  for (let frame = 0; frame < magSpec.length; frame++) {
    reconstructedSpec[frame] = [];
    for (let index = 0; index < magSpec[frame].length; index++) {
      let phase = math.complex(
        math.cos(phaSpec[frame][index][0]),
        math.sin(phaSpec[frame][index][0])
      );
      reconstructedSpec[frame].push(
        math.dotMultiply(
          math.dotPow(10, math.dotDivide(magSpec[frame][index][0], 20)),
          phase
        )
      );
    }
  }

  let reconstructedSpecTransposed = reconstructedSpec[0].map((_, colIndex) =>
    reconstructedSpec.map((row) => row[colIndex])
  );

  let signal = await griffinLim(reconstructedSpecTransposed, numIterations=100, windowLength=512, hopLength=256);
  callback(signal)

  return signal;
}
