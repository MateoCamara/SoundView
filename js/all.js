/* Global variables and objects */

// Login data
const CLIENT_ID = "J3QbU7A9Mt9wQbqYMRo9";
const CLIENT_SECRET = "TEWsO3ETlZ8aDPuvWYfqPyhYo97sl5COg9xEz4mO";
const REDIRECT_URL = "https://adrianbalda.github.io/Intellimixer/";
let AUTHORIZATION_CODE;
let accessToken;
let loginRedirected = false;
let loginSuccessful = 0;

// Tutorial
const tutorialPopup = document.getElementById('tutorial-popup');
const openTutorialPopup = document.getElementById('open-tutorial-popup');
const closeTutorialButton = document.getElementById('close-tutorial');
// Tutorial FIN

// Snackbar
let snackbar = document.getElementById("snackbar");
let snackbarMessage = document.getElementById("snackbarMessage");
// Snackbar FIN

//Menú Hamburguesa
const menuButton = document.querySelector(".menu-button");
const menuOptions = document.querySelector(".menu-options");
const menuContainer = document.querySelector(".menu-container");
const menuHamburguesa = document.getElementById("menuHamburguesa");
let menuVisible = false;
const queryForm = document.getElementById("query-form");
const uploadVAEs = document.getElementById('upload-vaes-div');
// Menú Hamburguesa FIN

// Transformacion de sonidos
const applyEffectsButton = document.getElementById('applyEffectsButton');
const gainControl = document.getElementById('gainControl');
const gainValue = document.getElementById('gainValue');
const speedControl = document.getElementById('speedControl');
const speedValue = document.getElementById('speedValue');
const effectsPopup = document.getElementById('effectsPopup');
// Transformacion de sonidos FIN

// Para el espectrograma
const heatmapColors = [
  '#000000', '#250066', '#4B0096', '#7200B6', '#9C00E6',
  '#C601E6', '#ED01C2', '#FF00A2', '#FF007F', '#FF005A',
  '#FF0036', '#FF0012', '#FF1100', '#FF2B00', '#FF4600',
  '#FF6000', '#FF7A00', '#FF9400', '#FFAF00', '#FFC900',
  '#FFE300', '#FFFD00', '#E8FF00', '#D3FF00', '#BEFF00',
  '#A9FF00', '#94FF00', '#7FFF00', '#6EFF00', '#5DFF00',
  '#4CFF00', '#3BFF00', '#2AFF00', '#1AFF00', '#09FF00',
  '#00FF0A', '#00FF1B', '#00FF2C', '#00FF3D', '#00FF4E',
  '#00FF5F', '#00FF70', '#00FF81', '#00FF92', '#00FFA3',
  '#00FFB4', '#00FFC5', '#00FFD6', '#00FFE7', '#00FFF8',
  '#00E6FF', '#00D2FF', '#00BDFF', '#00A8FF', '#0093FF',
  '#007EFF', '#0069FF', '#0055FF', '#0040FF', '#002BFF',
  '#0016FF', '#0000FF'
];

let dBData = [];
let maxValSoundsDB = [];
const canvasSpectrogram = document.getElementById('spectro');
const ctxSpectrogram = canvasSpectrogram.getContext('2d');
// Para el espectrograma

// Variational Autoencoder stuff
let encoderModel = undefined;
let decoderModel = undefined;
let model;
let sampleValuesInLatentSpace = [];
let sessionId = 0;
let numFiles = 3;
let latentSpaceLoaded = false;
let mu_latent_spaces = [];
let log_variance_latent_spaces = [];
let waveformImages = document.getElementById("waveform");
let idGeneratedAudio = -1;
let waveformAux;

// Audio stuff
let audio_manager = new AudioManager();
let MONO_MODE = true;

// Repeated strings
const DISPLAY_NONE = "none";
const DISPLAY_BLOCK = "block";
// Repeated strings

// Sounds and content
let default_query = "footstep";
let minDuration = 1;
let maxDuration = 2;
let sounds = [];
let selectedSound;
let currentSoundSamples = [];
let sampledSounds = [];
let extra_descriptors = undefined;
let map_features = undefined;
let n_pages = 3;
let n_pages_received = 0;
let all_loaded = false;
let last_selected_sound_id = undefined;
let soundsWaveforms = [];

// t-sne and xy map
let max_tsne_iterations = 500;
let current_it_number = 0;
let epsilon = 10;
let perplexity = 10;
let tsne = undefined;
let max_xy_iterations = 50;
let map_xy_x_max = undefined;
let map_xy_x_min = undefined;
let map_xy_y_max = undefined;
let map_xy_y_min = undefined;

// Canvas and display stuff
const popup = document.getElementById("loginPopup");
const overlay = document.querySelector(".overlay");
const soundInfoBox = document.getElementById('sound_info_box');
const switchButton = document.getElementById("switch-images-button");
let canvasWaveform = document.getElementById('waveform-generated');
let progressContainer = document.getElementById('progressContainer');
let canvasProgress = document.getElementById('canvasProgress');
let ctxProgress = canvasProgress.getContext('2d');
let ctxWaveform = canvasWaveform.getContext('2d');
let canvas = document.querySelector("canvas");
let ctx = canvas.getContext("2d");
canvasProgress.width = progressContainer.offsetWidth;
canvasProgress.height = 100;
let w = window.innerWidth;
let h = window.innerHeight;
let default_point_modulation = 0.6;
let disp_scale = Math.min(w, h);
let center_x = undefined; // Set in start()
let center_y = undefined; // Set in start()
let disp_x_offset;
let disp_y_offset;
let zoom_factor = undefined; // Set in start()
let rotation_degrees = undefined; // Set in start()
let min_zoom = 0.2;
let max_zoom = 10;
const PI = Math.PI;
const PI2 = 2 * Math.PI;

/* Setup and app flow functions */

function start() {
  // stop all audio
  audio_manager.stopAllBufferNodes();
  // get map descriptors
  setMapDescriptor();
  // update axis labels
  update_axis_labels();
  // Sounds
  sounds = [];
  sampledSounds = [];
  mu_latent_spaces = [];
  sampleValuesInLatentSpace = [];
  log_variance_latent_spaces = [];
  waveformImages.innerHTML = "";
  n_pages_received = 0;
  all_loaded = false;
  idGeneratedAudio = -1;

  // Canvas
  w = window.innerWidth;
  h = window.innerHeight;
  canvas.width = window.innerWidth;
  canvas.height = window.innerHeight;
  canvas.addEventListener("mousedown", onMouseDown, false);
  canvas.addEventListener("mouseup", onMouseUp, false);
  canvas.addEventListener("mouseout", onMouseOut, false);
  canvas.addEventListener("wheel", onWheel, false);
  center_x = 0.5;
  center_y = 0.5;
  zoom_factor = 1.0;
  rotation_degrees = 0;

  // Display stuff
  if (w >= h) {
    disp_x_offset = (w - h) / 2;
    disp_y_offset = 0.0;
  } else {
    disp_x_offset = 0.0;
    disp_y_offset = (h - w) / 2;
  }
  current_it_number = 0;

  dBData = [];
  maxValSoundsDB = [];
  soundsWaveforms = [];

  getSounds();
}

function getSounds(){
  //this is in online mode
  let query = document.getElementById("query_terms_input").value;

  // Search sounds in Freesound and start loading them
  if (query == undefined || query == "") {
    query = default_query;
  }
  let url =
  "https://freesound.org/apiv2/search/text/?query=" +
  query +
  "&group_by_pack=0" +
  "&filter=duration:[" +
  minDuration +
  "+TO+" +
  maxDuration +
  "]&page_size=" +
  numFiles +
  "&fields=id,previews,name,analysis,url,username,images";

  if(accessToken){
    url = url + "&token=" + accessToken.access_token + "&page=2";
  }

  loadJSON(
    function (data) {
      load_data_from_fs_json(data);
    },
    url,
    accessToken?.access_token
  );  
}

document.addEventListener('DOMContentLoaded', function() {
  closeTutorialButton.addEventListener('click', function() {
      tutorialPopup.style.display = 'none';
  });

  openTutorialPopup.addEventListener('click', function() {
    tutorialPopup.style.display = 'flex';
  });
});

// Efectos de Sonido
applyEffectsButton.addEventListener('click', () => {
  effectsPopup.style.display = DISPLAY_BLOCK;
  overlay.style.display = DISPLAY_BLOCK;
  speedControl.value = 1;
  speedValue.textContent = `x1.0`;
  gainControl.value = 1;
  gainValue.textContent = `1.0`;

  if(currentSoundSamples.length){
      let audioData = new Float32Array(currentSoundSamples);
      canvasProgress.width = progressContainer.offsetWidth;
      canvasProgress.height = window.innerWidth * 0.1;
      drawWaveform(audioData, canvasProgress, ctxProgress);    
  }
});

overlay.addEventListener('click', () => {
  effectsPopup.style.display = 'none';
  overlay.style.display = 'none';
  if (audioSource) {
      audioSource.stop(); // Si se esta reproduciendo el audio al salir de la transformacion de sonidos, se para.
  }    
});

effectsPopup.addEventListener('click', (event) => {
  event.stopPropagation();
});
// Efectos de Sonido FIN

switchButton.addEventListener("click", function() {
  canvasWaveform.style.display = (canvasWaveform.style.display === DISPLAY_NONE) ? DISPLAY_BLOCK : DISPLAY_NONE;
  canvasSpectrogram.style.display = (canvasSpectrogram.style.display === DISPLAY_NONE) ? DISPLAY_BLOCK : DISPLAY_NONE;

  if (canvasWaveform.style.display === DISPLAY_NONE) {
    switchButton.innerHTML = '<img src="https://img.icons8.com/external-smashingstocks-glyph-smashing-stocks/30/FFFFFF/external-bar-graph-shopping-and-commerce-smashingstocks-glyph-smashing-stocks.png" alt="external-bar-graph-shopping-and-commerce-smashingstocks-glyph-smashing-stocks"/>';
  } else {
    switchButton.innerHTML = '<img src="https://img.icons8.com/external-ayo-icons-royyan-wijaya/30/FFFFFF/external-waveform-audio-video-line-ayo-icons-royyan-wijaya.png" alt="external-waveform-audio-video-line-ayo-icons-royyan-wijaya"/>';
  }
});

window.addEventListener("load", async function () {
  AUTHORIZATION_CODE = getCodeFromURL();
  if(loginRedirected){
    popup.style.display = DISPLAY_NONE;
    overlay.style.display = DISPLAY_NONE;
    if(!accessToken){
      accessToken = await getAccessToken();
      getUserInfo(accessToken.access_token, function (userName) {
        showUser(userName);
      });
    }
    getSounds();
    loginRedirected = false;
  }
});

function showUser(userName) {
  const userNameElement = document.getElementById("userName");
  const username_container = document.getElementById("username-container");
  const login_container = document.getElementById("login-container");

  login_container.style.display = DISPLAY_NONE;

  username_container.style.display = DISPLAY_BLOCK;
  userNameElement.textContent = userName;
}

// Login popup
window.onload = function () {
  overlay.style.display = DISPLAY_BLOCK;

  popup.style.display = DISPLAY_BLOCK;
};

let formSubmitHandler = function formSubmitHandler(event) {
  event.preventDefault();
  start();
  hideForm();
  hideUploadVAEs()
};
queryForm.onsubmit = formSubmitHandler;

menuButton.addEventListener("click", function (event) {
  event.stopPropagation();
  if (menuVisible) {
    hideMenu();
  } else {
    showMenu();
  }
});

menuContainer.addEventListener("click", function (event) {
  event.stopPropagation();
});

document.addEventListener("click", function () {
  if (menuVisible) {
    hideMenu();
  }
});

document.getElementById('new-sound').addEventListener('click', function(event) {
  event.preventDefault();
  showForm();
  hideMenu();
});

document.getElementById('send-sound-request').addEventListener('click', formSubmitHandler);

document.getElementById('upload-vaes').addEventListener('click', function(event) {
  event.preventDefault();
  showUploadVAEs();
  hideMenu();
});

function showSnackbar(message) {
  snackbarMessage.innerHTML = message;
  snackbar.className = "show";
  setTimeout(function() { 
    if (snackbar.className === "show") {
      snackbar.className = snackbar.className.replace("show", ""); 
    }
  }, 10000);
}

function closeSnackbar() {
  snackbar.className = snackbar.className.replace("show", "");
}

function changeAxisAttribute() {
  all_loaded = false;
  current_it_number = 0;
  setMapDescriptor();
  update_axis_labels();
  center_x = 0.5;
  center_y = 0.5;
  zoom_factor = 1.0;
  rotation_degrees = 0;
}

window.requestAnimFrame = (function () {
  // This is called when code reaches this point
  return (
    window.requestAnimationFrame ||
    window.webkitRequestAnimationFrame ||
    window.mozRequestAnimationFrame ||
    function (callback) {
      window.setTimeout(callback, 1000 / 60);
    }
  );
})();

// Add the number of variables in latent space to html select tag items
function initLantentSpaceVariableSelector(latentSpaceDimension) {
  let xSelector = document.getElementById("x_axis_map_descriptors_selector");
  let ySelector = document.getElementById("y_axis_map_descriptors_selector");
  for (let i = 1; i <= latentSpaceDimension; i++) {
    let optX = document.createElement("option");
    optX.value = String(i);
    optX.innerHTML = "Dimension " + i;
    let optY = document.createElement("option");
    optY.value = String(i);
    optY.innerHTML = "Dimension " + i;
    xSelector.appendChild(optX);
    ySelector.appendChild(optY);
  }
  ySelector.getElementsByTagName("option")[1].selected = "selected"; // we use dim 1 and 2 in intiliazation
}

(async function init() {
  // This is called when code reaches this point
  window.addEventListener("keydown", onKeyDown, false);
  window.addEventListener("keyup", onKeyUp, false);
  // get encoder tensorflow model
  encoderModel = await tf.loadLayersModel(
    "https://mateocamara.com/intellimixer/models_tfg/encoder_model/model.json"
  );
  logInfo("Loaded default encoder.");
  // // get decoder tensorflow model
  decoderModel = await tf.loadLayersModel(
    "https://mateocamara.com/intellimixer/models_tfg/decoder_model/model.json"
  );
  logInfo("Loaded default decoder.");
  // Add the number of variables in latent space to html select tag items
  initLantentSpaceVariableSelector(encoderModel.outputShape[0][1]);
  setMapDescriptor(true);
  update_axis_labels();
})();

(function loop() {
  // This is called when code reaches this point
  // Get sound's xy position and scale it smoothly to create an animation effect
  if (mu_latent_spaces.length === numFiles && !all_loaded) {
    let x = mu_latent_spaces.map((ls) => ls[map_features[0]]);
    let y = mu_latent_spaces.map((ls) => ls[map_features[1]]);

    map_xy_x_max = Math.max.apply(null, x);
    map_xy_x_min = Math.min.apply(null, x);
    map_xy_y_max = Math.max.apply(null, y);
    map_xy_y_min = Math.min.apply(null, y);

    for (i in sounds) {
      sounds[i].computed_x =
        (x[i] - map_xy_x_min) / (map_xy_x_max - map_xy_x_min);
      sounds[i].computed_y =
        1 - (y[i] - map_xy_y_min) / (map_xy_y_max - map_xy_y_min);
    }

    for (i in sampledSounds) {
      sampledSounds[i].computed_x =
        (sampledSounds[i].latentSpace[map_features[0]] - map_xy_x_min) /
        (map_xy_x_max - map_xy_x_min);
      sampledSounds[i].computed_y =
        1 -
        (sampledSounds[i].latentSpace[map_features[1]] - map_xy_y_min) /
          (map_xy_y_max - map_xy_y_min);
    }

    all_loaded = true;
    latentSpaceLoaded = true;
    console.log("Loaded map with " + sounds.length + " sounds");
  }
  if (all_loaded && current_it_number <= max_xy_iterations) {
    for (i in sounds) {
      let sound = sounds[i];
      sound.x =
        sound.computed_x *
          Math.pow(100, current_it_number / max_xy_iterations - 1) +
        0.5 * (1 - Math.pow(100, current_it_number / max_xy_iterations - 1)); // Smooth position at the beginning
      sound.y =
        sound.computed_y *
          Math.pow(100, current_it_number / max_xy_iterations - 1) +
        0.5 * (1 - Math.pow(100, current_it_number / max_xy_iterations - 1)); // Smooth position at the beginning
    }

    for (i in sampledSounds) {
      let sound = sampledSounds[i];
      sound.x =
        sound.computed_x *
          Math.pow(100, current_it_number / max_xy_iterations - 1) +
        0.5 * (1 - Math.pow(100, current_it_number / max_xy_iterations - 1)); // Smooth position at the beginning
      sound.y =
        sound.computed_y *
          Math.pow(100, current_it_number / max_xy_iterations - 1) +
        0.5 * (1 - Math.pow(100, current_it_number / max_xy_iterations - 1)); // Smooth position at the beginning
    }
    current_it_number += 1;
  }
  if (all_loaded) {
    draw();
  }
  requestAnimFrame(loop);
})();

/* Sounds stuff */

function SoundFactory(id, preview_url, analysis, url, name, username, image) {
  this.x = Math.random();
  this.y = Math.random();
  this.rad = 15;
  this.mod_position = Math.random();
  this.mod_inc = 0.1;
  this.mod_amp = default_point_modulation;
  this.selected = false;

  this.id = id;
  this.preview_url = preview_url;
  this.analysis = analysis;

  let color = rgbToHex(255, 255, 255);
  this.rgba = color;

  this.url = url;
  this.name = name;
  this.username = username;
  this.image = image;
  this.generated = false;
}

function SoundFactoryGeneratedAudios(id, waveform, x, y, latentSpace) {
  this.id = id;
  this.waveform = waveform;
  this.x = x;
  this.y = y;
  this.rad = 15;
  this.mod_amp = default_point_modulation;
  this.mod_position = Math.random();
  this.mod_inc = 0.1;
  this.selected = false;
  this.generated = true;
  this.latentSpace = latentSpace;
}

function load_data_from_fs_json(data) {
  let max = data.results.length;
  let i = 0;
  try{
    if(!data.results.length){
      let querySoundsRequested = document.getElementById("query_terms_input").value;
      if(!querySoundsRequested){
        querySoundsRequested = default_query;
      }
      const message = `Couldn´t find "${querySoundsRequested}" sounds in Freesound.`;
      showSnackbar(message + " Please, try again with another sound!");
      throw new Error(message)
    }
    let interval = setInterval(function () {
      let sound_json = data.results[i];
      let sound = new SoundFactory(
          (id = sound_json.id),
          (preview_url =
            sound_json.audio || sound_json.previews["preview-hq-mp3"]),
          (analysis = sound_json.analysis),
          (url = sound_json.url),
          (name = sound_json.name),
          (username = sound_json.username),
          // (image = sound_json.images || sound_json.images.spectral_m)
          (image = [sound_json.images.spectral_m, sound_json.images.waveform_m])
      );
      sounds.push(sound);

      getWaveformFromPreview(function (waveform) {
        let adjustedWaveform = adjustAudioToExpectedSize(waveform, 22050);
        sessionId += 1;
        soundsWaveforms.push(waveform);

        // TODO
        const signal = tf.tensor1d(adjustedWaveform);
        const frameLength = 512;
        const frameStep = Math.floor(frameLength / 4); // Si no, con 248 se obtiene 64x256
        const fftLength = 256;
        const signalTransformed = tf.signal
          .stft(signal, frameLength, frameStep, fftLength)
          .arraySync();
        const finalData = [];
        signalTransformed.forEach((element) => {
          const data = element.slice(0, 256);
          finalData.push(data);
        });
        
        if (finalData.length > 64) {
          finalData.splice(64);
        } else if (finalData.length < 64) {
          const numRowsToAdd = 64 - finalData.length;
          const lastRow = finalData[finalData.length - 1];
          for (let i = 0; i < numRowsToAdd; i++) {
            finalData.push(lastRow);
          }
        }
        let mcltspecTransposed = finalData[0].map((_, colIndex) =>
          finalData.map((row) => row[colIndex])
        );
        
        // Para el espectrograma
        const stftAbs = tf.abs(mcltspecTransposed);

        // Convertir a decibelios (escala logarítmica)
        const stftDb = tf.tidy(() => {
          const minAmp = tf.max(stftAbs).div(1e6).clipByValue(1, Infinity);
          const logSpec = tf.log(stftAbs.add(minAmp));
          return logSpec.mul(10).div(tf.log(tf.scalar(10)));
        });
                    
        // Obtener los valores de amplitud del espectrograma
        const stftData = stftDb.arraySync();
        dBData.push(stftData);
        maxValSoundsDB.push(tf.max(stftData).dataSync()[0]);
        // Para el espectrograma

        let mclt2Dspec = spectrogram(mcltspecTransposed);
        let { mu_latent_space, log_variance_latent_space } =
          encodeAudio(mclt2Dspec);
        let latent_space_size = mu_latent_space.length;
        mu_latent_spaces.push(mu_latent_space);
        log_variance_latent_spaces.push(log_variance_latent_space);
      }, sound.preview_url);
      i++;
      if (i == max) {
        clearInterval(interval);
      }
    }, 1000);
  }catch(error){
    console.log(error.message)
  }
}

function checkSelectSound(x, y) {
  applyEffectsButton.style.display = DISPLAY_BLOCK;
  let min_dist = 9999;
  let selected_sound = false;
  let distancesArray = [];
  let spectro_selected_sound = [];
  let max_value_spectro = [];
  let waveform_selected_Sound = [];
  currentSoundSamples = [];
  for (i in sounds) {
    let sound = sounds[i];
    let dist = computeEuclideanDistance(sound.x, sound.y, x, y);
    if (dist < min_dist) {
      min_dist = dist;
      selected_sound = sound;
      selectedSound = sound;
      spectro_selected_sound = dBData[i];
      max_value_spectro = maxValSoundsDB[i];
      waveform_selected_Sound = soundsWaveforms[i];
      currentSoundSamples = soundsWaveforms[i];
    }
    distancesArray.push(dist);
  }
  for (i in sampledSounds) {
    let sound = sampledSounds[i];
    let dist = computeEuclideanDistance(sound.x, sound.y, x, y);
    if (dist < min_dist) {
      min_dist = dist;
      selected_sound = sound;
      selectedSound = sound;
    }
  }
  if (min_dist < 0.02) {
    selectSound(selected_sound, spectro_selected_sound, max_value_spectro, waveform_selected_Sound);
  } else {
    let dim1LatentSpace = x * (map_xy_x_max - map_xy_x_min) + map_xy_x_min;
    let dim2LatentSpace =
      (y - 1) * -(map_xy_y_max - map_xy_y_min) + map_xy_y_min;

    let quadraticDistancesArray = distancesArray.map((dist) => dist ** 2);
    let minDistancePointIndex = distancesArray.indexOf(
      Math.min(...distancesArray)
    );

    let multilateralDim = [];
    let aux = [];
    for (let i = 0; i < mu_latent_spaces.length; i++) {
      if (i === minDistancePointIndex) continue;
      for (j in mu_latent_spaces[i]) {
        aux.push(
          (quadraticDistancesArray[i] -
            quadraticDistancesArray[minDistancePointIndex] -
            mu_latent_spaces[i][j] ** 2 +
            mu_latent_spaces[minDistancePointIndex][j] ** 2) /
            (2 *
              (mu_latent_spaces[minDistancePointIndex][j] -
                mu_latent_spaces[i][j]))
        );
      }
      multilateralDim.push(aux);
      aux = [];
    }

    multilateralDim = multilateralDim[0];

    // var result_array = [];
    // for (i = 0; i < 20; i++) {
    //     result_array.push((mu_latent_spaces[0][i] + mu_latent_spaces[1][i] + mu_latent_spaces[2][i]) / 3)
    // }
    // multilateralDim = result_array;

    multilateralDim[map_features[0]] = dim1LatentSpace;
    multilateralDim[map_features[1]] = dim2LatentSpace;

    //  multilateralDim = [-0.70845157,  0.37552756, -1.4685907 ,  0.18429044,  0.12212666,
    //  1.820008  ,  0.5608268 , -0.78984535,  0.5061305 ,  0.07133903,
    //  0.35714608, -3.453202  , -0.7022973 ,  0.9746381 ,  1.4159285 ,
    // -1.7573189 ,  0.7665267 ,  0.4188148 ,  1.1196271 , -0.40932375];

    let latent_space_size = mu_latent_spaces[0].length;
    let decoder_tensor = tf.tensor(multilateralDim);
    let predicted_spectogram_tensor = decoderModel.predict(
      tf.reshape(decoder_tensor, (shape = [1, latent_space_size]))
    );
    let predicted_spectrogram = predicted_spectogram_tensor.arraySync();
    predicted_spectrogram = predicted_spectrogram[0];
    convertPredictedSpectrogramIntoAudio(
      function (audio) {
        idGeneratedAudio += 1;
        let sound = new SoundFactoryGeneratedAudios(
          (id = idGeneratedAudio),
          (waveform = audio),
          (x = x),
          (y = y),
          (latentSpace = multilateralDim)
        );
        sampledSounds.push(sound);
        playGeneratedSound(sound.waveform);
        showGeneratedSoundInfo(sound.waveform);
      },
      predicted_spectrogram,
      "2D",
      512
    );
  }
}

function selectSound(selected_sound, spectro_selected_sound, max_value_spectro, waveform_selected_Sound) {
  selected_sound.selected = true;
  selected_sound.mod_amp = 5.0;
  if (MONO_MODE) {
    audio_manager.stopAllBufferNodes();
  }
  if (selected_sound.generated) {
    playGeneratedSound(selected_sound.waveform);
    showGeneratedSoundInfo(selected_sound.waveform);
  } else {
    audio_manager.loadSound(selected_sound.id, selected_sound.preview_url);
    showSoundInfo(selected_sound, spectro_selected_sound, max_value_spectro, waveform_selected_Sound);
  }
  last_selected_sound_id = selected_sound.id;
  selected_sound.selected = false;
  selected_sound.mod_amp = default_point_modulation;
}

function finishPlayingSound(sound_id) {
  let sound = getSoundFromId(sound_id);
  sound.selected = false;
  sound.mod_amp = default_point_modulation;
}

function selectSoundFromId(sound_id) {
  let sound = getSoundFromId(sound_id);
  selectSound(sound);
}

function getSoundFromId(sound_id) {
  for (i in sounds) {
    let sound = sounds[i];
    if (sound.id == parseInt(sound_id)) {
      return sound;
    }
  }
}

function showSoundInfo(sound, spectro_selected_sound, max_value_spectro, waveform_selected_Sound) {
  switchButton.innerHTML = '<img src="https://img.icons8.com/external-ayo-icons-royyan-wijaya/30/FFFFFF/external-waveform-audio-video-line-ayo-icons-royyan-wijaya.png" alt="external-waveform-audio-video-line-ayo-icons-royyan-wijaya"/>';
  switchButton.style.display = DISPLAY_BLOCK
  canvasSpectrogram.style.display = DISPLAY_NONE

  let html = "";
  html +=
    sound.name +
    ' by <a href="' +
    sound.url +
    '" target="_blank">' +
    sound.username +
    "</a>";
    const soundInfoContent = document.getElementById("sound_info_content");

  // Para el espectrograma. Representar el espectrograma en el lienzo (canvas)
  for (let y = 0; y < spectro_selected_sound.length; y++) {
    for (let x = 0; x < spectro_selected_sound[y].length; x++) {
      const colorValue = spectro_selected_sound[y][x];
      const color = getColorFromValue(colorValue, max_value_spectro);
      ctxSpectrogram.fillStyle = color;
      ctxSpectrogram.fillRect(x, y, 1, 1);
    }
  }
// Para el espectrograma
  let audioData = new Float32Array(waveform_selected_Sound);
  canvasWaveform.style.display = DISPLAY_BLOCK;
  canvasWaveform.width = soundInfoBox.offsetWidth;
  canvasWaveform.height = 100;
  soundInfoBox.appendChild(canvasWaveform);
  drawWaveform(audioData, canvasWaveform, ctxWaveform);
  soundInfoContent.innerHTML = html;

  //Para la transformacion
  canvasProgress.width = progressContainer.offsetWidth;
  canvasProgress.height = 100;
  drawWaveform(audioData, canvasProgress, ctxProgress);
}

function drawWaveform(data, waveCanvas, waveCtx) {
  waveCtx.clearRect(0, 0, waveCanvas.width, waveCanvas.height);
  waveCtx.strokeStyle = 'yellow';
  waveCtx.beginPath();
  waveCtx.moveTo(0, (1 + data[0]) * waveCanvas.height / 2);
  for (let i = 1; i < data.length; i++) {
    waveCtx.lineTo(i * waveCanvas.width / data.length, (1 + data[i]) * waveCanvas.height / 2);
  }
  waveCtx.stroke();
}

function updateProgress(progress) {
  const progressColor = 'orange';

  drawWaveform(currentSoundSamples, canvasProgress, ctxProgress)

  // Pintar progreso en color naranja
  ctxProgress.strokeStyle = progressColor;
  ctxProgress.beginPath();
  ctxProgress.moveTo(0, (1 + currentSoundSamples[0]) * canvasProgress.height / 2);

  const progressIndex = Math.floor(progress * currentSoundSamples.length);
  for (let i = 1; i <= progressIndex; i++) {
    const x = i * canvasProgress.width / currentSoundSamples.length;
    const y = (1 + currentSoundSamples[i]) * canvasProgress.height / 2;
    ctxProgress.lineTo(x, y);
  }
  ctxProgress.stroke();
}

// Para el espectrograma
function getColorFromValue(value, maxValue) {
  const normalizedValue = value / maxValue;
  const colorIndex = Math.floor(normalizedValue * (heatmapColors.length - 1));
  return heatmapColors[colorIndex];
}

function showGeneratedSoundInfo(waveform) {
  waveformAux = waveform;
  let html = "";
  html += `<button onclick="downloadAudio(waveformAux)">Download generated audio!</button>`;
  document.getElementById("sound_info_box").innerHTML = html;
}

function downloadAudio(waveform) {
  const waveformArray = new Float32Array(waveform);
  const wavBytes = getWavBytes(waveformArray.buffer, {
    isFloat: true, // floating point or 16-bit integer
    numChannels: 1,
    sampleRate: 22050,
  });
  const wav = new Blob([wavBytes], { type: "audio/wav" });
  const downloadLink = document.createElement("a");
  downloadLink.href = URL.createObjectURL(wav);
  downloadLink.setAttribute("download", "my-audio.wav");
  downloadLink.click();
}

function playGeneratedSound(waveform) {
  let wave = waveform.map((x) => x * 20);
  let myArrayBuffer = context.createBuffer(1, 22050 * 3, 22050);
  let wf = new Float32Array(wave);
  let nowBuffering = myArrayBuffer.getChannelData(0);
  for (var i = 0; i < myArrayBuffer.length; i++) {
    if (i >= wf.length) {
      nowBuffering[i] = 0;
    } else {
      nowBuffering[i] = wf[i];
    }
  }
  var source = context.createBufferSource();
  source.buffer = myArrayBuffer;
  source.connect(context.gainNode);
  source.start();
}

function setMapDescriptor(init) {
  let x_descriptor =
    document.getElementById("x_axis_map_descriptors_selector").value - 1;
  let y_descriptor =
    document.getElementById("y_axis_map_descriptors_selector").value - 1;
  if (init) {
    y_descriptor = 1;
  }
  map_features = [x_descriptor, y_descriptor];
}

/* Drawing */

function draw() {
  ctx.clearRect(0, 0, w, h);
  ctx.globalCompositeOperation = "lighter";
  for (i in sounds) {
    let sound = sounds[i];
    let disp_x, disp_y;
    [disp_x, disp_y] = normCoordsToDisplayCoords(sound.x, sound.y);

    let z = i;
    z == sounds.length - 1 ? (z = 0) : (z = 1 + parseInt(i));
    let nextSound = sounds[z];
    let nextX, nextY;
    [nextX, nextY] = normCoordsToDisplayCoords(nextSound.x, nextSound.y);

    ctx.fillStyle = "#ffffff";
    ctx.strokeStyle = "#ffffff";
    ctx.beginPath(); // Start a new path
    ctx.moveTo(disp_x, disp_y); // Move the pen to (30, 50)
    ctx.lineTo(nextX, nextY); // Draw a line to (150, 100)
    ctx.stroke();
    ctx.closePath();

    if (!sound.selected) {
      ctx.fillStyle = sound.rgba;
      ctx.strokeStyle = sound.rgba;
    } else {
      ctx.fillStyle = "#ffffff";
      ctx.strokeStyle = "#ffffff";
    }

    if (last_selected_sound_id == sound["id"]) {
      ctx.fillStyle = "#0000ff";
      ctx.strokeStyle = "#0000ff";
    }

    ctx.beginPath();
    ctx.arc(
      disp_x,
      disp_y,
      sound.rad * zoom_factor * Math.pow(0.9, zoom_factor),
      0,
      Math.PI * 2,
      true
    );
    ctx.fill();
    ctx.closePath();

    ctx.beginPath();
    ctx.arc(
      disp_x,
      disp_y,
      (sound.rad + 5 + sound.mod_amp * Math.cos(sound.mod_position)) *
        zoom_factor *
        Math.pow(0.9, zoom_factor),
      0,
      Math.PI * 2,
      true
    );
    ctx.stroke();
    ctx.closePath();

    sound.mod_position += sound.mod_inc;
  }

  for (i in sampledSounds) {
    let sound = sampledSounds[i];
    let disp_x, disp_y;
    [disp_x, disp_y] = normCoordsToDisplayCoords(sound.x, sound.y);

    ctx.fillStyle = "#ffff00";
    ctx.strokeStyle = "#ffff00";

    ctx.beginPath();
    ctx.arc(
      disp_x,
      disp_y,
      sound.rad * zoom_factor * Math.pow(0.9, zoom_factor),
      0,
      Math.PI * 2,
      true
    );
    ctx.fill();
    ctx.closePath();

    ctx.beginPath();
    ctx.arc(
      disp_x,
      disp_y,
      (sound.rad + 5 + sound.mod_amp * Math.cos(sound.mod_position)) *
        zoom_factor *
        Math.pow(0.9, zoom_factor),
      0,
      Math.PI * 2,
      true
    );
    ctx.stroke();
    ctx.closePath();

    sound.mod_position += sound.mod_inc;
  }
}

// form submit event handler
(function () {
  let formSubmitHandler = function formSubmitHandler(event) {
    event.preventDefault();
    start();
  };
  document.getElementById("query-form").onsubmit = formSubmitHandler;
})();

// axis text label drawing
function update_axis_labels() {
  // update the text boxes
  let feature1 = 1 + +map_features[0];
  let feature2 = 1 + +map_features[1];
  document.getElementById("x_axis_box").innerHTML = "Dimension " + feature1;
  document.getElementById("y_axis_box").innerHTML = "Dimension " + feature2;
}
