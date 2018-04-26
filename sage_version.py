import re
import urllib2

# Obtain the different Sage versions
def get_all_version_names(mirror_url, idx = None, distribution = '{{cookiecutter.travis_ubuntu_version}}'):
    if idx is None:
        idx = 0
    else:
        idx = int(idx)
    site = urllib2.urlopen(mirror_url).read()
    ans = re.findall('(sage-([0-9]*(?:\.[0-9]*)*)-%s.tar.bz2)'%distribution, site)
    all_version_names = []
    for fname, ver in ans:
        if fname not in all_version_names:
            all_version_names.append(fname)
    return all_version_names[idx]
