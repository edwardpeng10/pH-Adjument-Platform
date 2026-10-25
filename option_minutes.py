"""v16: robot minutes of the titrations (no pools) for each start-prompt combination, on product B (the 12 Oct
product), with tips 5-8 on and off - the numbers the two start prompts quote.    python tools/option_minutes.py"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import scen16 as scen
if __name__ == '__main__':
    for eight in (1, 0):
        row = {}
        for name, neut, base in (('acid only', 0, 0), ('acid + VI', 1, 0), ('acid + base', 0, 1), ('all three', 1, 1)):
            s = scen.run('B', 12, answers=dict(neut=neut, base_on=base, pool_lanes='0'), override={'Use_8_tips': eight})
            blocks = ('acid_c1', 'acid_c2', 'val', 'tris_val')
            row[name] = sum(s.block_sim.get(b, 0) for b in blocks) / 60 + (0 if not s.errors else float('nan'))
        print(f"tips {'1-8' if eight else '1-4'}: " + ', '.join(f'{k} {v:.1f} min' for k, v in row.items())
              + f" | VI adds {row['acid + VI'] - row['acid only']:.1f}, base adds {row['acid + base'] - row['acid only']:.1f}")
