import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex
import Image, ImageDraw, ImageFont

cells = ['Gm12878', 'H1hesc', 'Helas3', 'Hepg2', 'Huvec', 'K562', 'Nhek']
colors = [rgb2hex(plt.cm.Set1(c)) for c in np.linspace(0,1,9)]

im = Image.new('RGBA', (120, (len(colors)+2)*20), '#FFFFFF')

draw = ImageDraw.Draw(im)
font = ImageFont.truetype('/usr/share/fonts/truetype/freefont/FreeSans.ttf', 14)

s = 0
scx = 20
scy = 20
for label, color in zip(cells, colors):
	draw.rectangle([(scx, scy+s), (scx+20, scy+s+20)], fill=color)
	draw.text((scx+30, scy+3+s), label, fill='#000000', font=font)
	s += 20

im.save('legend.png', 'PNG')


