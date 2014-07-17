#!/bin/csh -f

set file = $1

setenv IMAGIC_BATCH 1
echo "! "
echo "! "
echo "! ====================== "
echo "! IMAGIC ACCUMULATE FILE "
echo "! ====================== "
echo "! "
echo "! "
echo "! IMAGIC program: em2em ------------------------------------------------"
echo "! "
/programs/x/imagic/leschziner/stand/em2em.e <<EOF
TIF
SPI
SINGLE_FILE
NO
NO
${file:r}
${file:r}.spi
spi
LINUX
EOF

