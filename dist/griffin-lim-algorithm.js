function istft_(realpart, impart) {
    const nFrames = 64 //realPart.length;
    let parArray = []
    let imparArray = []

    for (let n = 0; n < nFrames; n++) {
        let waveform = tf.spectral.irfft(tf.complex(realpart[n], impart[n])) //tf.real(tf.ifft(tf.complexSlice(X, [n, 0], [1, N])));
        let waveform_array = waveform.arraySync()[0]

        if (n % 2 === 0) {
            parArray = parArray.concat(waveform_array)
        } else {
            imparArray = imparArray.concat(waveform_array)
        }

    }

    parArray = parArray.concat(new Array(256).fill(0).flat())
    imparArray.splice(0, 0, new Array(256).fill(0))
    imparArray = imparArray.flat(1)

    let result = math.add(parArray, imparArray)

    return result;
}

function angle(x) {
    const real = tf.real(x);
    const imag = tf.imag(x);
    return tf.atan2(imag, real);
}


async function griffinLim(magnitudeSpectrogram, numIterations=100, windowLength=512, hopLength=256) {
    let phase = tf.randomNormal(tf.tensor(magnitudeSpectrogram).shape);
    phase = phase.arraySync()
    let result;
    let coseno
    let seno
    let complexo

    for (let i = 0; i < numIterations; i++) {
        coseno = math.cos(phase)
        seno = math.sin(phase)
        complexo = tf.complex(coseno, seno)
        result = istft_(coseno, seno)
        if (i < numIterations - 1) {
            let stftResult = tf.signal.stft(tf.tensor(result), 512, 256);

            let phase = angle(stftResult).arraySync();
        }
    }

    return tf.tensor(result).arraySync();
}