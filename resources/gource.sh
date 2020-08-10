#!/bin/bash

# Don't forget to use Ctrl - + to speed up the video (otherwise it takes forever)
gource --title "APFEL++" --hide-filenames --key -a 1 --path ../. --highlight-all-users -o - | ffmpeg -y -r 60 -f image2pipe -vcodec ppm -i - -vcodec libx264 -preset ultrafast -pix_fmt yuv420p -crf 20 -threads 0 -bf 0 progress.mp4
