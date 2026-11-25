// Record the README demo of the three.js replay (viz/replay/index.html, built by `python -m viz.build`).
//   node viz/record_demo.mjs   -> out/viz/demo/*.webm + replay.png
// then: ffmpeg -ss 4 -t 11 -i out/viz/demo/*.webm -c:v libx264 -pix_fmt yuv420p -crf 23 -r 25 docs/media/replay.mp4
//       (GIF: fps=12, scale=960, 96-colour palette from the MP4)
// Headless Chromium with SwiftShader WebGL; fails if the page logs a console error. ALL NUMBERS ARE MODELLED.
import { createRequire } from 'node:module';
import { mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';
const require = createRequire(import.meta.url);
const { chromium } = require('playwright');

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const out = path.join(root, 'out', 'viz', 'demo');
mkdirSync(out, { recursive: true });
const b = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
const ctx = await b.newContext({ viewport: { width: 1280, height: 720 }, recordVideo: { dir: out, size: { width: 1280, height: 720 } } });
const page = await ctx.newPage();
const errors = [];
page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
await page.goto('file://' + path.join(root, 'viz', 'replay', 'index.html'));
await page.waitForTimeout(16000);
await page.screenshot({ path: path.join(out, 'replay.png') });
await ctx.close();
await b.close();
if (errors.length) { console.error(errors); process.exit(1); }
console.log('recorded to', out);
