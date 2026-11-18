# Procedural nuclear explosion - volumetric ray marcher in C

A single C file that renders one frame of a nuclear blast (flash -> fireball -> shock-raised base surge -> rising cap that rolls into a vortex ring -> dust stem) as a 3D volume. It is not a fluid simulation. Each frame is built from time-dependent shapes deformed by noise, then ray-marched with self-shadowing and temperature-based emission.

## Build

```bash
gcc -O3 -march=native -ffast-math -fopenmp nuke.c -o nuke -lm
```

## Run

```bash
./nuke <t> <out.raw> [W H CX GY SC ROLL]
```

| Argument | Meaning | Default |
|---|---|---|
| `t` | time in seconds; the sequence spans 0–4 s | — |
| `out.raw` | output file | — |
| `W H` | image size in px | 758 630 |
| `CX GY` | screen position of ground zero (centre x, ground line y) | 379 557.5 |
| `SC` | pixels per world unit (zoom) | 200 |
| `ROLL` | camera roll in radians | 0 |

The output is raw `float32 [H][W][4]`: linear HDR RGB (premultiplied, values can exceed 1) plus alpha (cloud opacity). stdout prints the fireball's screen x, y and the cap radius in px.

For camera shake, vary `CX`, `GY` and `ROLL` from frame to frame.

Frame `i` at 30 fps: `./nuke $(echo "$i/30" | bc -l) frame_$i.raw`

## Algorithm

1. **Timeline (`setpar`).** Closed-form curves over `t` set the fireball radius `0.56(1−e^(−t/0.1))`, the rise height, the cap's ellipsoid axes, the sphere→ring morph `m = smoothstep(0.75, 2.5, t)`, the stem's height and width, the base surge front `0.25 + 1.28(1−e^(−(t−0.05)/0.5))`, and the core and surface temperatures.
2. **Shapes as signed distance functions.**
   - **Cap:** an ellipsoid smoothly merged with a torus, with the underside smoothly hollowed out.
   - **Stem:** a wobbling column that flares at the bottom and top.
   - **Base surge:** a rolling elliptical wall plus an inner layer, with its radius and height varied by noise around the ring.
3. **Chaotic surfaces.** The distance field is displaced by inverted Worley noise at two scales (the cauliflower lumps), billow noise (Σ|Perlin|, fine puffs) and a time-varying Perlin domain warp (boiling). The noise coordinates follow the flow: they roll around the torus core, scale with growth, rise and swirl in the stem, and expand radially in the base. Density is `clamp(−sdf/edge, 0, 1)`.
4. **Temperature -> emission.** Temperature depends on depth inside the cloud, height, and bulge vs. crevice, with a soot mask and flicker. It maps to colour through a black-body-style HDR ramp.
5. **Lighting.**
   - A 120×80×120 grid of transmittance toward the key light (Beer–Lambert).
   - A one-sample self-shadow probe for each billow.
   - The fireball as a point light on the dust, with `1/(1+1.3d²)` falloff.
   - Dim red ambient light.
6. **Ray marching (OpenMP, one ray per pixel).** An orthographic, slightly tilted camera. Rays skip empty space using a conservative bound, then take fixed jittered steps of 0.0085. Compositing is front to back: `α = 1−e^(−σΔ)`, `C += Tr·α·(albedo·light + emission)`, `Tr *= 1−α`. Rays exit early when `Tr < 0.004`. Rays that hit the ground plane pick up glow from the fireball.
