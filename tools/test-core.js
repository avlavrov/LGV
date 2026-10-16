// Проверка ядра на демо-данных: число ридов, MD-восстановление референса против FASTA, известные варианты.
'use strict';
const fs = require('fs'), path = require('path');
const C = require('../src/core.js');
const D = path.join(__dirname, '..', 'demo');
const truth = JSON.parse(fs.readFileSync(path.join(D, 'truth.json')));
const blob = f => new Blob([fs.readFileSync(path.join(D, f))]);
let fails = 0;
const ok = (cond, msg) => { console.log((cond ? 'OK   ' : 'FAIL ') + msg); if (!cond) fails++; };

(async () => {
  const bam = await new C.BamFile(new C.FileSource(blob('demo.bam')), new C.FileSource(blob('demo.bam.bai'))).open();
  ok(bam.refs.length === 2 && bam.refs[0].name === 'demo' && bam.refs[0].len === 40000, 'заголовок: ' + JSON.stringify(bam.refs) + ' образцы ' + bam.samples);
  const fa = await new C.IndexedFasta(new C.FileSource(blob('demo.fa')), new C.FileSource(blob('demo.fa.fai'))).open();

  let t = Date.now();
  const all = await bam.fetch(0, 0, 40000);
  const all2 = await bam.fetch(1, 0, 2000);
  ok(all.length + all2.length === truth.reads, `все риды: ${all.length}+${all2.length} = ${truth.reads} (${Date.now() - t} мс)`);
  ok(new Set(all.map(r => r.name + (r.flag & 64))).size === all.length, 'нет дублей при объединении чанков');

  const at = await bam.fetch(0, 12000, 12001);
  ok(at.length === truth.overlap_12000_12001, `перекрывают 12000: ${at.length} = ${truth.overlap_12000_12001}`);
  const at2 = await bam.fetch(1, 1000, 1200);
  ok(at2.length === truth.overlap_demo2_1000_1200, `demo2:1000-1200: ${at2.length} = ${truth.overlap_demo2_1000_1200}`);

  // MD → референс совпадает с FASTA
  const ws = 0, refArr = new Uint8Array(40000);
  for (const r of all) C.fillRefFromMD(r, refArr, ws);
  const faSeq = await fa.fetch('demo', 0, 40000);
  let filled = 0, bad = 0; for (let i = 0; i < 40000; i++) if (refArr[i]) { filled++; if (refArr[i] !== faSeq[i]) bad++; }
  ok(bad === 0 && filled > 30000, `референс из MD: заполнено ${filled}, расхождений с FASTA ${bad}`);
  const fa2 = await fa.fetch('demo2', 995, 1005);
  ok(fa2.length === 10, 'FASTA второго контига: ' + String.fromCharCode(...fa2));

  // варианты
  const [sp, sref, salt] = truth.variants.hetSnv;
  const cov = C.computeCoverage(all.filter(r => !(r.flag & 1024)), 0, 40000);
  const idx = { A: 0, C: 1, G: 2, T: 3 };
  const c12 = cov.slice(sp * 6, sp * 6 + 6);
  ok(c12[idx[salt]] > 5 && c12[idx[sref]] > 5, `гетерозиготный SNV ${sref}>${salt} @${sp}: ${Array.from(c12)}`);
  const [hp, , halt] = truth.variants.homSnv; const c15 = cov.slice(hp * 6, hp * 6 + 6);
  ok(c15[idx[halt]] / (c15[0] + c15[1] + c15[2] + c15[3]) > 0.9, `гомозиготный SNV @${hp}: ${Array.from(c15)}`);
  ok(cov[18001 * 6 + 5] > 5, 'делеция 3 п.н.: ' + cov[18001 * 6 + 5]);
  const ins = all.filter(r => C.cigarString(r).includes('3I')).length;
  ok(ins > 5, 'риды со вставкой 3 п.н.: ' + ins);
  const refAt = p => faSeq[p];
  const mism = all.filter(r => r.pos <= sp && r.end > sp).flatMap(r => C.computeFeatures(r, refAt)).filter(f => f.t === 0 && f.pos === sp).length;
  ok(mism > 5, 'несовпадения в позиции SNV: ' + mism);
  const mid = cov[26500 * 6] + cov[26500 * 6 + 1] + cov[26500 * 6 + 2] + cov[26500 * 6 + 3];
  ok(mid > 5 && mid < 25, 'покрытие внутри крупной гетерозиготной делеции ~половина: ' + mid);
  const orient = {}; for (const r of all) { const o = C.pairOrientation(r); orient[o] = (orient[o] || 0) + 1; }
  console.log('ориентации пар', orient, 'границы вставки', C.insertSizeBounds(all));
  ok(orient.other >= 24, 'химерные пары');
  const rows = C.packRows(all.slice(0, 3000), 2);
  let overlap = 0; for (const row of rows) for (let i = 1; i < row.length; i++) if (row[i].pos < row[i - 1].end + 2) overlap++;
  ok(overlap === 0, `раскладка: ${rows.length} рядов, наложений ${overlap}`);
  console.log(fails ? `\n${fails} ПРОВАЛ(ОВ)` : '\nВсе проверки пройдены');
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(1); });
