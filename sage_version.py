import re
import urllib2

# Obtain the different Sage versions
def get_all_version_names(mirror_url, idx = None, distribution = 'Ubuntu_18.04-x86_64'):
    if idx is None:
        idx = 0
    else:
        idx = int(idx)
    all_version_names = []
    for subdir in ["", "old/"]:
        site = urllib2.urlopen(mirror_url + subdir).read()
        ans = re.findall('(sage-([0-9]*(?:\.[0-9]*)*)-%s.tar.bz2)'%distribution, site)
        for fname, ver in ans:
            if fname not in all_version_names:
                all_version_names.append(subdir + fname)
    return all_version_names[idx]
