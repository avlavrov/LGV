import sys, json
from playwright.sync_api import sync_playwright
url = 'file:///home/claude/cgv/dist/cgv.html'
out = sys.argv[1] if len(sys.argv) > 1 else '/tmp/claude-0/-home-claude/05d33790-9911-556a-819e-597227a9ba65/scratchpad'
steps = json.loads(sys.argv[2]) if len(sys.argv) > 2 else []
with sync_playwright() as p:
    b = p.chromium.launch(args=['--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--ignore-gpu-blocklist'])
    for scheme in ['light', 'dark']:
        pg = b.new_page(viewport={'width': 1280, 'height': 760}, color_scheme=scheme)
        logs = []
        pg.on('console', lambda m: logs.append(m.type + ': ' + m.text))
        pg.on('pageerror', lambda e: logs.append('PAGEERROR: ' + str(e)))
        pg.goto(url)
        pg.wait_for_timeout(1500)
        for s in steps:
            if s[0] == 'locus':
                pg.fill('#locus', s[1]); pg.press('#locus', 'Enter'); pg.wait_for_timeout(900)
            elif s[0] == 'js':
                pg.evaluate(s[1]); pg.wait_for_timeout(600)
            elif s[0] == 'click':
                pg.mouse.click(s[1], s[2]); pg.wait_for_timeout(400)
            elif s[0] == 'sel':
                pg.click(s[1]); pg.wait_for_timeout(400)
        print(scheme, 'status:', pg.inner_text('#status'))
        print(scheme, 'locus:', pg.input_value('#locus'))
        for l in logs: print('  ', l)
        pg.screenshot(path=f'{out}/cgv-{scheme}.png')
        if scheme == 'light':
            pg2 = b.new_page(viewport={'width': 400, 'height': 760}, color_scheme=scheme)
            pg2.goto(url); pg2.wait_for_timeout(1500)
            pg2.screenshot(path=f'{out}/cgv-phone.png')
    b.close()
