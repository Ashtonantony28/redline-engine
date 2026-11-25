// Headless renderer for cad/render.py: node render.mjs <workdir> <outdir> <width> <height> <scale> <vendor_dir>
// Serves <workdir> at https://render.local/ and three.js from <vendor_dir> (same URLs as the CDN), renders every
// view of scene.json and writes <outdir>/<view name>.png. Prints one JSON line {ok, errors, views}.
import { createRequire } from 'module';
import fs from 'fs';
import path from 'path';
const require = createRequire(import.meta.url);
const { chromium } = require('playwright');

const [work, outdir, w, h, scale, vendor] = process.argv.slice(2);
const CDN = 'https://cdn.jsdelivr.net/npm/';
const report = { ok: false, errors: [], views: [] };
const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader',
  '--ignore-gpu-blocklist'] });
try {
  const page = await browser.newPage({ viewport: { width: +w, height: +h }, deviceScaleFactor: +scale });
  page.on('console', (m) => { if (m.type() === 'error') report.errors.push(m.text()); });
  page.on('pageerror', (e) => report.errors.push(String(e)));
  await page.route('**/*', async (route) => {
    const url = route.request().url();
    let local = null, type = 'application/octet-stream';
    if (url.startsWith('https://render.local/')) local = path.join(work, decodeURIComponent(url.slice(21)));
    else if (url.startsWith(CDN)) local = path.join(vendor, url.slice(CDN.length));
    if (local && fs.existsSync(local)) {
      if (local.endsWith('.js')) type = 'application/javascript';
      else if (local.endsWith('.html')) type = 'text/html';
      else if (local.endsWith('.json')) type = 'application/json';
      else if (local.endsWith('.glb')) type = 'model/gltf-binary';
      return route.fulfill({ status: 200, contentType: type, body: fs.readFileSync(local) });
    }
    report.errors.push('blocked ' + url);
    return route.abort();
  });
  await page.goto('https://render.local/index.html');
  await page.waitForFunction(() => window.__ready === true, null, { timeout: 120000 });
  const scene = JSON.parse(fs.readFileSync(path.join(work, 'scene.json'), 'utf8'));
  for (let i = 0; i < scene.views.length; i++) {
    const info = await page.evaluate(async (k) => await window.renderView(k), i);
    const file = path.join(outdir, scene.views[i].name + '.png');
    await page.locator('canvas').screenshot({ path: file, timeout: 900000 });  // software GL: large scenes are slow
    report.views.push({ name: scene.views[i].name, file, ...info });
  }
  report.ok = report.errors.length === 0;
} catch (e) {
  report.errors.push('render failed: ' + e);
} finally {
  await browser.close();
}
console.log(JSON.stringify(report));
process.exit(report.ok ? 0 : 1);
