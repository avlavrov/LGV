// Builds one self-contained HTML: template + i18n + core + app + demo data (base64).
'use strict';
const fs = require('fs'), path = require('path');
const R = p => path.join(__dirname, '..', p);
const b64 = f => fs.readFileSync(R('demo/' + f)).toString('base64');
const demo = `const DEMO = ${JSON.stringify({ bam: b64('demo.bam'), bai: b64('demo.bam.bai'), fa: b64('demo.fa'), fai: b64('demo.fa.fai') })};`;
const safe = s => s.replace(/<\/script/gi, '<\\/script');
let html = fs.readFileSync(R('src/index.html'), 'utf8');
const BUILD_ID = new Date().toISOString(), BUILD_DATE = BUILD_ID.slice(0, 10);
const MANE_VER = (fs.readFileSync(R('data/mane.grch38.tsv'), 'utf8').match(/^#(?:CGV|LGV) MANE \S+ (v[\d.]+)/m) || [, ''])[1];
html = html.replace('/*DEMO*/', () => demo).replace('/*I18N*/', () => safe(fs.readFileSync(R('src/i18n.js'), 'utf8'))).replace('/*CORE*/', () => safe(fs.readFileSync(R('src/core.js'), 'utf8') + '\n' + fs.readFileSync(R('src/compress.js'), 'utf8') + '\n' + fs.readFileSync(R('src/htscodecs.js'), 'utf8') + '\n' + fs.readFileSync(R('src/cram.js'), 'utf8')))  .replace('/*APP*/', () => safe((fs.readFileSync(R('src/genes.js'), 'utf8') + '\n' + fs.readFileSync(R('src/tracks.js'), 'utf8') + '\n' + fs.readFileSync(R('src/signal.js'), 'utf8') + '\n' + fs.readFileSync(R('src/app.js'), 'utf8') + '\n' + fs.readFileSync(R('src/signalview.js'), 'utf8')).replace('__BUILD_DATE__', BUILD_DATE).replace('__BUILD_ID__', BUILD_ID).replace('__MANE__', MANE_VER)))
  // MANE transcripts (data/mane.grch38.tsv), packed (tools/pack-genes.js) + base64: read by src/genes.js at start
  .replace('/*CYTO*/', () => safe(fs.readFileSync(R('data/cytoband.grch38.tsv'), 'utf8')))
  .replace('/*MANE*/', () => require('./pack-genes.js').packGenes(fs.readFileSync(R('data/mane.grch38.tsv'), 'utf8')).toString('base64'));
// every inline script must compile: a syntax error would leave the page blank
{
  const vm = require('vm');
  for (const m of html.matchAll(/<script(?![^>]*type="text\/plain")[^>]*>([\s\S]*?)<\/script>/g)) {
    try { new vm.Script(m[1]); } catch (e) { console.error('Build stopped: a script does not compile:', e.message); process.exit(1); }
  }
}
fs.mkdirSync(R('dist'), { recursive: true });
// for claude.ai: no <html>/<head> of its own — the platform adds them
fs.writeFileSync(R('dist/artifact.html'), html);
// for local use and any web server: a full document with a strict CSP
const local = '<!doctype html>\n<html lang="ru">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
  + '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; script-src \'unsafe-inline\'; style-src \'unsafe-inline\'; img-src data: blob:; connect-src https: http://localhost:* http://127.0.0.1:*; worker-src blob:">\n'
  + '<style>body{margin:0}</style>\n</head>\n<body>\n' + html + '\n</body>\n</html>\n';
fs.writeFileSync(R('dist/LGV.html'), local);
// what the update check of saved copies compares (published next to LGV.html)
const VERS = fs.existsSync(R('data/versions.json')) ? JSON.parse(fs.readFileSync(R('data/versions.json'), 'utf8')) : {};
fs.writeFileSync(R('dist/version.json'), JSON.stringify({ build: BUILD_ID, mane: MANE_VER, gencode: VERS.gencode || '' }) + '\n');
console.log('dist/LGV.html', (local.length / 1024).toFixed(0), 'KB; dist/artifact.html', (html.length / 1024).toFixed(0), 'KB');
