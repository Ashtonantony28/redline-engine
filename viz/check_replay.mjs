// Headless check of the Redline replay page (called by `python -m viz.build --check` and the M10 tests).
// usage: node check_replay.mjs <html> <screenshot.png> <settle_ms> <gap_ms> <width> <height> <cdn> <vendor_dir>
// CDN requests are answered from <vendor_dir> (files fetched by viz/build.py with curl, same URLs and versions)
// so the check does not depend on the browser reaching the CDN; anything else off-page is aborted and reported.
// Prints one JSON line: {ok, console_errors, page_errors, theta0, theta1, frames0, frames1, ...}.
import { createRequire } from 'module';
import fs from 'fs';
import path from 'path';
const require = createRequire(import.meta.url);
const { chromium } = require('playwright');

const [html, shot, settle, gap, w, h, cdn, vendor] = process.argv.slice(2);
const report = { ok: false, console_errors: [], page_errors: [], blocked: [], served_from_vendor: 0 };
const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'] });
try {
  const page = await browser.newPage({ viewport: { width: +w, height: +h } });
  page.on('console', (m) => { if (m.type() === 'error') report.console_errors.push(m.text()); });
  page.on('pageerror', (e) => report.page_errors.push(String(e)));
  await page.route('**/*', async (route) => {
    const url = route.request().url();
    if (url.startsWith('file:') || url.startsWith('data:')) return route.continue();
    if (url.startsWith(cdn + '/')) {
      const local = path.join(vendor, url.slice(cdn.length + 1));
      if (fs.existsSync(local)) {
        report.served_from_vendor++;
        return route.fulfill({ status: 200, contentType: 'application/javascript', body: fs.readFileSync(local) });
      }
    }
    report.blocked.push(url);
    return route.abort();
  });
  await page.goto('file://' + path.resolve(html));
  await page.waitForSelector('body[data-ready="1"]', { timeout: 90000 });
  await page.waitForTimeout(+settle);
  const s0 = await page.evaluate(() => ({ theta: window.__replay.theta, frames: window.__replay.frames }));
  await page.waitForTimeout(+gap);
  const s1 = await page.evaluate(() => ({ theta: window.__replay.theta, frames: window.__replay.frames }));
  report.theta0 = s0.theta; report.theta1 = s1.theta; report.frames0 = s0.frames; report.frames1 = s1.frames;
  report.webgl = await page.evaluate(() => { const c = document.querySelector('#view canvas'); return !!(c && (c.getContext('webgl2') || c.getContext('webgl'))); });
  await page.screenshot({ path: shot });
  report.screenshot = shot;
  report.ok = report.console_errors.length === 0 && report.page_errors.length === 0 && report.blocked.length === 0 &&
              s1.frames > s0.frames && s1.theta !== s0.theta;
} catch (e) {
  report.page_errors.push('check failed: ' + e);
} finally {
  await browser.close();
}
console.log(JSON.stringify(report));
process.exit(report.ok ? 0 : 1);
