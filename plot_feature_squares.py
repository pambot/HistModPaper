import sys
import numpy as np
import cPickle as pickle
import Image, ImageDraw


keys = ['H', 'B', 'C']

def get_zeroed_indexes(fs_fname):
    hist_scores = pickle.load(open('predict4/results/' + fs_fname, 'rb'))
    feat = hist_scores.keys()
    fs_scores = np.array(hist_scores.values())
    ind = dict.fromkeys(keys, None)
    ind['H'] = np.array([i for i,v in enumerate(feat) if v[0]=='H' and v[-1]=='3'])
    ind['B'] = np.array([i for i,v in enumerate(feat) if v[0]=='H' and v[-1]=='1'])
    ind['C'] = np.array([i for i,v in enumerate(feat) if v[0]=='E'])
    return {k:fs_scores[ind[k]] for k in keys}


def draw_circles(k, (cx, cy), color='k', w=20, s=5):
    cx1, cy1 = cx, cy
    hs, vs = s, s
    for i, ind in enumerate(fs_dict[k]):
        try:
            val = fs_dict[k][i]
            if val:
                draw.ellipse([(cx, cy), (cx+w, cy+w)], fill=color)
            else:
                draw.ellipse([(cx, cy), (cx+w, cy+w)], fill='#d8dcd6')
        except IndexError:
            break
        cx += hs + w
        if (i + 1) % 10 == 0:
            cy += vs + w
            cx = cx1
    return


fs_dict = get_zeroed_indexes('histScoresCodingGm12878.pkl')

w, s = 20, 5
x_start = 10
im = Image.new('RGBA', ((w + s)*11*3, 600), '#FFFFFF')

draw = ImageDraw.Draw(im)

colors = [
'#738595',
'#03719c', 
'#c65102',
]

for k, color in zip(keys, colors):
    draw_circles(k, (x_start, 10), color=color, w=w, s=s)
    x_start += (w + s)*10 + w

im.save('figures/feature_select_example.png')



