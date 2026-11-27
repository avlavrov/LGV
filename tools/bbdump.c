// libBigWig bigBed dumper for tools/test-bigbed.js (independent check of tools/bigbed.js): prints every record.
// Build: gcc -O2 -DNOCURL -o bbdump tools/bbdump.c LIBBIGWIG/{bwRead,bwStats,bwValues,bwWrite,io}.c -ILIBBIGWIG -lz -lm
#include "bigWig.h"
#include <stdio.h>
#include <stdlib.h>
// usage: bbdump file.bb  -> prints every record: chrom start end rest  (all contigs, whole length)
int main(int argc, char **argv) {
  if (bwInit(1<<17) != 0) return 2;
  bigWigFile_t *bb = bbOpen(argv[1], NULL);
  if (!bb) { fprintf(stderr, "open failed\n"); return 1; }
  char *sql = bbGetSQL(bb); fprintf(stderr, "fields=%d defined=%d sql=%s\n", bb->hdr->fieldCount, bb->hdr->definedFieldCount, sql ? "yes" : "no"); free(sql);
  for (int i = 0; i < bb->cl->nKeys; i++) {
    bbOverlappingEntries_t *o = bbGetOverlappingEntries(bb, bb->cl->chrom[i], 0, bb->cl->len[i], 1);
    if (!o) { fprintf(stderr, "query failed %s\n", bb->cl->chrom[i]); return 1; }
    for (uint32_t j = 0; j < o->l; j++) printf("%s\t%u\t%u\t%s\n", bb->cl->chrom[i], o->start[j], o->end[j], o->str[j]);
    bbDestroyOverlappingEntries(o);
  }
  bwClose(bb); bwCleanup(); return 0;
}
