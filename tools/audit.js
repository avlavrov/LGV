// Security audit of CGV's own code. Run: node tools/audit.js (part of every check before a release).
// Fails on patterns that are never allowed, and lists the few places that need a human look so a
// reviewer sees every one of them. See SECURITY.md.
'use strict';
const fs = require('fs'), path = require('path');
const R = p => path.join(__dirname, '..', p);
const files = ['src/core.js', 'src/compress.js', 'src/htscodecs.js', 'src/cram.js', 'src/genes.js', 'src/tracks.js', 'src/app.js', 'src/i18n.js', 'src/index.html'];
let fails = 0;
const fail = (f, n, line, why) => { fails++; console.log(`FAIL ${f}:${n}  ${why}\n       ${line.trim().slice(0, 150)}`); };

// never allowed anywhere in our code
const FORBIDDEN = [
  [/\beval\s*\(/, 'eval runs text as code'],
  [/\bnew\s+Function\s*\(|[^.\w]Function\s*\(/, 'Function() runs text as code'],
  [/\bset(Timeout|Interval)\s*\(\s*['"`]/, 'a timer with a string runs text as code'],
  [/document\.write/, 'document.write'],
  [/\b(outerHTML|insertAdjacentHTML)\b/, 'HTML from strings outside the reviewed places'],
  [/\bimportScripts\s*\(|\bimport\s*\(/, 'loading code at run time'],
  [/\bWebSocket\b|\bEventSource\b|sendBeacon|XMLHttpRequest/, 'a network channel CGV does not use'],
  [/<script[^>]*\bsrc=|<link[^>]*\bhref=["']https?:/i, 'a remote script or stylesheet'],
  [/window\.(opener|parent|top)\b/, 'talking to other windows'],
  [/document\.cookie/, 'cookies'],
];
// allowed only in the listed files; every other use fails
const PLACES = [
  // a call of the browser's global fetch (methods named fetch on our own classes are not network calls)
  [code => /(^|[^.\w$])fetch\s*\(/.test(code.replace(/\basync\s+fetch\s*\(/g, '')),
    (f, code) => f === 'src/core.js' && /\bfetch\s*\(\s*this\.url/.test(code), 'network request outside UrlSource (only links the user gave are read)'],
  [code => /localStorage/.test(code), (f, code) => f === 'src/app.js' && /PROFILE_KEY/.test(code), 'browser storage other than the display profile'],
  [code => /https?:\/\/[^\s'"`)…]/.test(code),
    (f, code, line) => (f === 'src/app.js' && /^\s*(hg38|hg19|chm13):/.test(line)) || (f === 'src/i18n.js') || (f === 'src/index.html' && /placeholder=/.test(line)),
    'a fixed web address (allowed: the public reference presets only)'],
];
const reviewed = []; // innerHTML lines: each must pass every value through esc()

for (const f of files) {
  const lines = fs.readFileSync(R(f), 'utf8').split('\n');
  lines.forEach((line, i) => {
    const n = i + 1, code = line.replace(/(^|[^:\\])\/\/.*$/, '$1').replace(/\/\*.*?\*\//g, '');
    if (/^\s*(\/\/|\*)/.test(line) || !code.trim()) return;
    for (const [re, why] of FORBIDDEN) if (re.test(code)) fail(f, n, line, why);
    for (const [hit, allowed, why] of PLACES) if (hit(code) && !allowed(f, code, line)) fail(f, n, line, why);
    if (/innerHTML\s*=/.test(code)) {
      // every ${...} on the line is esc(...), a number/boolean we built, or an already-escaped fragment
      const parts = [...code.matchAll(/\$\{([^}]*)\}/g)].map(m => m[1].trim());
      const loose = parts.filter(x => !/^esc\(/.test(x) && !/^(cols\.length|rows|code|i|kind|on|act|val|\(?on \?|P\.cardPair|!P\.cardPair|cls \?|P\.cardPair \?)/.test(x) && !/^col\.title \?|^goMate \?|^linesHtml\(/.test(x));
      reviewed.push(`${f}:${n}`);
      if (loose.length) fail(f, n, line, 'innerHTML with a value not passed through esc(): ' + loose.join(', '));
    }
  });
}
// file readers report damage only as CgvError (a clear message for the user), never a bare Error
for (const f of ['src/core.js', 'src/compress.js', 'src/htscodecs.js', 'src/cram.js', 'src/tracks.js']) fs.readFileSync(R(f), 'utf8').split('\n').forEach((line, i) => {
  if (/throw\s+(new\s+(?!CgvError\b)\w*Error\b|['"`])/.test(line.replace(/\/\/.*$/, ''))) fail(f, i + 1, line, 'a file reader throws something other than CgvError');
});
// every reader of file data has a fuzzing campaign in the tests
{
  const tests = ['tools/test-fuzz.js', 'tools/test-cram.js'].map(f => fs.readFileSync(R(f), 'utf8')).join('\n');
  const need = { 'BAM records': /BAM records, damage/, 'BAM header': /BAM header, damage/, BAI: /BAI index, damaged/, BGZF: /BGZF block headers/, '.fai': /FASTA index \(\.fai\), damaged/, '.gzi': /FASTA \.gzi, damaged/, 'FASTA indexing': /FASTA indexed in the browser, damaged/, CRAM: /fuzzing records/, 'CRAM codecs': /fuzzing/, locus: /locus field/, profile: /pasted profile/, 'annotation tracks': /Annotation tracks \(BED/ };
  const miss = Object.keys(need).filter(k => !need[k].test(tests));
  if (miss.length) { fails++; console.log('FAIL no fuzzing campaign for: ' + miss.join(', ')); } else console.log('OK   every file reader has a fuzzing campaign: ' + Object.keys(need).join(', '));
}
// top-level names must not repeat across the scripts that share one page (core + cram + app + i18n)
{
  const seen = new Map();
  for (const f of ['src/i18n.js', 'src/core.js', 'src/compress.js', 'src/htscodecs.js', 'src/cram.js', 'src/genes.js', 'src/tracks.js', 'src/app.js']) for (const m of fs.readFileSync(R(f), 'utf8').matchAll(/^(?:const|let|var|class|function|async function) ([A-Za-z_$][\w$]*)/gm)) {
    if (seen.has(m[1])) { fails++; console.log(`FAIL ${f}: top-level name ${m[1]} is also declared in ${seen.get(m[1])}`); } else seen.set(m[1], f);
  }
  console.log('OK   top-level names checked across the page\'s scripts: ' + seen.size);
}
// the local build must carry the strict Content Security Policy
const dist = R('dist/cgv.html');
if (fs.existsSync(dist)) {
  const h = fs.readFileSync(dist, 'utf8');
  const csp = (h.match(/Content-Security-Policy" content="([^"]+)"/) || [])[1] || '';
  if (!/default-src 'none'/.test(csp) || /unsafe-eval/.test(csp) || !/script-src 'unsafe-inline'/.test(csp)) { fails++; console.log('FAIL dist/cgv.html: the Content Security Policy is missing or weakened: ' + csp); }
  else console.log('OK   dist/cgv.html carries the strict Content Security Policy');
}
// no third-party packages
if (fs.existsSync(R('package.json'))) { const pj = JSON.parse(fs.readFileSync(R('package.json'))); if (Object.keys(pj.dependencies || {}).length) { fails++; console.log('FAIL package.json lists runtime dependencies'); } }
console.log(`OK   innerHTML places, each checked for esc(): ${reviewed.join(', ')}`);
console.log(fails ? `\n${fails} problem(s)` : '\nAudit passed');
process.exit(fails ? 1 : 0);
