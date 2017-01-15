BEGIN {FS = "\t"}
{
  printf("%s%c", $1, "\t"); 
  for (i = 2; i <= NF; i += 4) printf ("%s%c", $i, i + 4 <= NF ? "\t" : "\n");
}

