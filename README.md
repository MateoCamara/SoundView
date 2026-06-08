# SoundView (Intellimixer)

An interactive **sound map** in the browser. Sounds are laid out in 2D by the
dimensions of a **Variational Autoencoder (VAE)** latent space — so similar-sounding
clips sit close together — and you can **generate brand-new sounds** by clicking on
empty regions of the map: the click position is turned into a latent vector,
decoded back into a spectrogram, and inverted to audio.

Runs entirely in the browser with TensorFlow.js — no server, no build step required.

## Demo (offline by default)

Just serve the folder and open it:

```bash
# any static server works; for example:
python -m http.server 8000
# then open http://localhost:8000
```

On load it reads a small bundled dataset (`files/latents.json`, 50 sounds with their
pre-computed VAE latents) and the local VAE models in `models/`, so it works with **no
login and no network calls** for the map itself.

- **Click a point** → listen to that sound (streamed from Freesound's public preview CDN).
- **Click an empty area** → generate a new sound at that latent position.
- **Click + drag** to pan, **scroll** to zoom, and use the **axis selectors** to map
  different latent dimensions to X/Y.

## How it works

- `models/encoder_model` / `models/decoder_model` — the VAE (latent dim 4), loaded with
  `tf.loadLayersModel`.
- `files/latents.json` — each sound's metadata plus its encoder output (`mu`). Built
  offline by running every clip through the encoder once, so the map loads instantly and
  needs no audio fetching at runtime (Freesound previews are CORS-restricted, so they
  can't be decoded client-side; they are only *played* via an `<audio>` element).
- Generation: latent vector → `decoderModel.predict` → spectrogram (magnitude + phase) →
  a single inverse STFT → waveform.

To rebuild `files/latents.json` (e.g. for more sounds), download the previews referenced
in `files/BlogPostDemo.json` locally and encode them with the same pipeline the app uses.

## Optional: live Freesound search

The bundled demo does not require Freesound. To enable live search, define a config
object **before** the app scripts load and never commit a real client secret (use the
OAuth2 implicit/PKCE flow so no secret is needed client-side):

```html
<script>
  window.FREESOUND_CONFIG = {
    clientId: "your-client-id",
    redirectUrl: "https://your-site/",
    enableLogin: true,
  };
</script>
```

## Credits

This project builds on **[Intellimixer](https://github.com/adrianBalda/Intellimixer) by
Adrián Balda** (the original sound-map concept, canvas UI and Freesound integration).
The VAE-based placement and the generate-a-sound pipeline were added by Mateo Cámara.
Sounds are from [Freesound](https://freesound.org); audio analysis/synthesis uses
[TensorFlow.js](https://www.tensorflow.org/js).

## License

See [LICENSE](LICENSE). The MIT license covers Mateo Cámara's modifications and additions;
the upstream Intellimixer it derives from is published without a license, so those portions
remain under their original author's copyright.
