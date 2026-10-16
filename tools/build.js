// Builds one self-contained HTML: template + i18n + core + app + demo data (base64).
'use strict';
const fs = require('fs'), path = require('path');
const R = p => path.join(__dirname, '..', p);
const b64 = f => fs.readFileSync(R('demo/' + f)).toString('base64');
const demo = `const DEMO = ${JSON.stringify({ bam: b64('demo.bam'), bai: b64('demo.bam.bai'), fa: b64('demo.fa'), fai: b64('demo.fa.fai') })};`;
const safe = s => s.replace(/<\/script/gi, '<\\/script');
let html = fs.readFileSync(R('src/index.html'), 'utf8');
html = html.replace('/*DEMO*/', () => demo).replace('/*I18N*/', () => safe(fs.readFileSync(R('src/i18n.js'), 'utf8'))).replace('/*CORE*/', () => safe(fs.readFileSync(R('src/core.js'), 'utf8'))).replace('/*APP*/', () => safe(fs.readFileSync(R('src/app.js'), 'utf8')));
fs.mkdirSync(R('dist'), { recursive: true });
// for claude.ai: no <html>/<head> of its own — the platform adds them
fs.writeFileSync(R('dist/artifact.html'), html);
// for local use and any web server: a full document with a strict CSP
const local = '<!doctype html>\n<html lang="ru">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
  + '<meta http-equiv="Content-Security-Policy" content="default-src \'none\'; script-src \'unsafe-inline\'; style-src \'unsafe-inline\'; img-src data: blob:; connect-src https: http://localhost:* http://127.0.0.1:*; worker-src blob:">\n'
  + '<style>body{margin:0}</style>\n</head>\n<body>\n' + html + '\n</body>\n</html>\n';
fs.writeFileSync(R('dist/cgv.html'), local);
console.log('dist/cgv.html', (local.length / 1024).toFixed(0), 'KB; dist/artifact.html', (html.length / 1024).toFixed(0), 'KB');
