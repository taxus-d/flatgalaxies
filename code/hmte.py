# -*- coding: utf-8 -*-
"""
    code.hmte
    ~~~~~~~~~

    Home-made template engine (HMTE)

    expand loop-like templates in sql queries
    see .tsql files in queries directory for examples

    :copyright: (c) 2019 by taxus-d.
    :license: MIT, see LICENSE for more details.
"""

import string
import re

import pandas as pd  # required by another part of project anyway
import numpy as np
import yaml


class FormatDict(dict):
    """
    An extension of python .format that preserves the patterns
    it wasn't told about in processed string
    """
    def __missing__(self, key):
        return "{" + key + "}"


_formatter = string.Formatter()


_statregexes = {
    'loop start': re.compile(r"(\s*)--\s*for\s*(\w+)\s*in\s*([^-]*)--\s*"),
    'end'       : re.compile(r"--\s*end\s*--")
}


def _trim_body_sniff_comma(body):
    """
    detect and trim comma at the end of sting if it is present
    """
    trimmer = re.compile(r"(,?)\s+$")
    has_comma = False
    if "," in re.findall(trimmer, body):
        has_comma = True
    body_trimmed = re.sub(trimmer, "", body)
    return body_trimmed, has_comma


def handle_loops(s, stats, sb, se, level):
    """
    Expand all loops
    The main routine in hmte currently

    May be generalized to arbitary structure handler in distant future
    """
    if len(stats) > 0:
        bottomstats = stats.loc[stats['level'] == level]
        begins = bottomstats.iloc[range(0, len(bottomstats), 2)]
        ends = bottomstats.iloc[range(1, len(bottomstats), 2)]
        try:
            assert len(begins) == len(ends)
            assert np.all(ends.type == 'end')
        except AssertionError:
            raise SyntaxError(
                """your `begins' dont match your `ends', check template throughly
                """)

        rep_body = ""
        cpos = sb

        for i in range(len(begins)):
            b, e = begins.iloc[i], ends.iloc[i]

            grs = _statregexes['loop start'].fullmatch(s[b.beg:b.end])

            indent = grs[1]
            loop_variable = grs[2]
            lis = yaml.load(grs[3], yaml.SafeLoader)

            idx_start = begins.index[i]
            idx_stop  = ends.index[i]

            # recurse !
            # (loops can be arbitary nested)
            body = handle_loops(s, stats.loc[idx_start+1:idx_stop-1],
                                b.end, e.beg, level+1)

            rep_body += s[cpos:b.beg]

            for i in range(len(lis)):
                v = lis[i]
                body_trimmed, has_comma_at_end = _trim_body_sniff_comma(body)
                optcomma = "," if (i < len(lis)-1 or level > 0) else ""
                optcaret = indent  # if i > 0 else ""
                # matching formatting is always tricky...

                mapping = FormatDict(**{loop_variable: v})
                rep_body += optcaret + _formatter.vformat(body_trimmed, (), mapping) + optcomma
            cpos = e.end
        rep_body += s[cpos:se]
    else:
        rep_body = s[sb:se]
    return rep_body


def find_loops(s, stats):
    loopstarts = re.finditer(_statregexes['loop start'], s)
    newrows = []
    for stat in loopstarts:
        newrows.append({
            "beg": stat.span()[0],
            "end": stat.span()[1],
            "type": "loop start"
        })
    if newrows == []:
        return stats
    else:
        return stats.append(newrows, ignore_index=True)


def find_ends(s, stats):
    ends = re.finditer(_statregexes['end'], s)
    newrows = []
    for stat in ends:
        newrows.append({
            "beg": stat.span()[0],
            "end": stat.span()[1],
            "type": "end"
        })
    if newrows == []:
        return stats
    else:
        return stats.append(newrows, ignore_index=True)


def assign_levels(stats):
    stats.level = 0
    currentlevel = 0
    for i in range(len(stats)):
        t = stats.loc[i]["type"]
        if t in ("loop start"):
            stats.loc[i, "level"] = currentlevel
            currentlevel += 1
        elif t in ("end"):
            currentlevel -= 1
            stats.loc[i, "level"] = currentlevel
    return stats


def expand_templates(t, **kw):
    mapping = FormatDict(**kw)
    query_ws = _formatter.vformat(t, (), mapping)

    s = query_ws
    stats = pd.DataFrame(columns=["beg", "end", "type", "level"])

    stats = find_loops(s, stats)
    stats = find_ends(s, stats)
    stats = stats.sort_values(by="beg")
    stats.index = range(len(stats))
    stats = assign_levels(stats)
    query = handle_loops(s, stats, 0, len(s), 0)
    return query
