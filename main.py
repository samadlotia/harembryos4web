import re
import os
from os.path import join as filejoin, isdir, isfile
import sys
import csv
import shutil

HAR_ID_RE = re.compile(r'(?:2x)?HAR\.?(\d+)')

ACTIVITY_DOMAIN_RE = re.compile(r'(.+)\s+\((\d+)\)')
def parse_activity_domains(s):
  if not s or s == 'none':
    return []
  pieces = s.split(';')
  domains = []
  for piece in pieces:
    piece = piece.strip()
    match = ACTIVITY_DOMAIN_RE.match(piece)
    if not match: continue
    name = match.group(1)
    num_imgs = int(match.group(2))
    domains.append((name, num_imgs))
  return domains

def parse_har_id(s):
  match = HAR_ID_RE.match(s)
  if not match:
    return None
  else:
    return int(match.group(1))

GENE_AND_LOCATION_RE = re.compile(r'^(\S+)\s*\((\-?\d+)\)$')
def parse_genes(s):
  if not s or s == 'None': return None
  pieces = s.split(',')
  genes = []
  for piece in pieces:
    piece = piece.strip()
    match = GENE_AND_LOCATION_RE.match(piece)
    if not match:
      continue
    gene = match.group(1)
    location = int(match.group(2))
    genes.append((gene, location))
  return genes

def parse_bracketed_genes(s, har_name):
  genes = parse_genes(s)
  if not genes: return None
  upstreams = filter(lambda g: g[1] < 0, genes)
  if not upstreams:
    print 'Warning: no upstream genes found for', har_name
    return None
  downstreams = filter(lambda g: g[1] > 0, genes)
  if not downstreams:
    print 'Warning: no downstream genes found for', har_name
    return None
  nearest_upstream = max(upstreams, key = lambda g: g[1])
  nearest_downstream = min(downstreams, key = lambda g: g[1])
  return [nearest_upstream[0], nearest_downstream[0]]

def parse_aliases(s):
  if not s:
    return []
  pieces = s.split(';')
  return map(lambda s: parse_har_id(s.strip()), pieces)

def parse_har_csv(csv_path):
  hars = {}
  input_file = open(csv_path, 'r')
  input = csv.DictReader(input_file)
  for row in input:
    if not row['ID']:
      continue
    har_id = parse_har_id(row['ID'])
    if not har_id:
      print "Warning: skipping unrecognized HAR ID in '%s': %s" % (csv_path, row['ID'])
      continue
    if not har_id in hars:
      hars[har_id] = {
        'name'               : row['ID'],
        'species'            : {},
        'aliases'            : parse_aliases(row['Aliases']),
        'species-difference' : row['Human-Chimp Difference'],
      }
    if not hars[har_id].get('bracketed-genes'):
      hars[har_id]['bracketed-genes'] = parse_bracketed_genes(row['Genes within 1 Mb (distance to TSS)'], row['ID'])

    hars[har_id]['species'][row['Species'].lower()] = {
      'genome-coords'               : row['Coordinates (hg19 or panTro4)'],
      'consistent-activity-domains' : parse_activity_domains(row['Consistent Activity Domains (# pos)']),
      'suggestive-activity-domains' : parse_activity_domains(row['Suggestive Activity Domains (# pos)']),
      'expression'                  : row['Expression'],
      'stage'                       : row['Stage'],
      'imgs'                        : {},
    }
  input_file.close()
  return hars

ALLOWED_IMG_RE = re.compile(r'(?P<har>\d+)_(?P<species>[^\W\d_]{2})\d+_(?P<num>\d+)L\.tif$')

IMG_SPECIES_TO_CSV = {
  'hg': 'human',
  'pt': 'chimp'
}

def find_imgs(hars, img_dir_path, skipped):
  for dirname, _, filenames in os.walk(img_dir_path):
    for filename in filenames:
      filepath = filejoin(dirname, filename)
      match = ALLOWED_IMG_RE.match(filename)
      if not match:
        continue
      har = int(match.group('har'))
      if not har in hars:
        if skipped: print "Skipping -- unknown HAR '%d': %s" % (har, filepath)
        continue
      img_species = match.group('species')
      species = IMG_SPECIES_TO_CSV.get(img_species)
      if not species:
        if skipped: print "Skipping -- unknown species '%s': %s" % (species, filepath)
        continue
      img_num = int(match.group('num'))
      hars[har]['species'][species]['imgs'][img_num] = filepath

OUTPUT_DIR_PREFIX = 'website'
def make_output_dir_name():
  dir_name = OUTPUT_DIR_PREFIX
  n = 1
  while isdir(dir_name):
    n += 1
    dir_name = OUTPUT_DIR_PREFIX + ' ' + str(n)
  return dir_name

THUMBNAIL_WIDTH = 280
def calc_thumbnail_size(img_size):
  w, h = img_size
  tw = THUMBNAIL_WIDTH 
  th = tw * h / w
  return (tw, th)

def convert_img(original_path, target_prefix):
  import Image # Python caches import statements
  try:
    img = Image.open(original_path)
  except Exception, e:
    print "Error: could not open image because: %s\n  %s" % (str(e), original_path)
    return
  img.save(target_prefix + '.tif', format = 'TIFF') # copy full tiff image
  thumbnail_img = img.resize(calc_thumbnail_size(img.size), Image.ANTIALIAS)
  thumbnail_img.save(target_prefix + '-small.png', format = 'PNG')

def convert_imgs(hars, output_dir):
  img_dir = filejoin(output_dir, 'imgs')
  os.mkdir(img_dir)
  for har in hars:
    har_dir = filejoin(img_dir, str(har))
    if not isdir(har_dir): os.mkdir(har_dir)
    for species in hars[har]['species']:
      species_dir = filejoin(har_dir, species)
      if not isdir(species_dir): os.mkdir(species_dir)
      for img_num, original_path in hars[har]['species'][species]['imgs'].iteritems():
        target_path = filejoin(species_dir, str(img_num))
        convert_img(original_path, target_path)

GENOME_URL_PREFIXES = {
  'human': 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=',
  'chimp': 'http://genome.ucsc.edu/cgi-bin/hgTracks?org=Chimp&db=panTro4&position='
}

def setup_web_files(hars, output_dir):
  from jinja2 import Environment, FileSystemLoader
  env = Environment(loader = FileSystemLoader('.'))
  template_input = env.get_template('index-template.html')
  output_file = open(filejoin(output_dir, 'index.html'), 'w')
  template_input.stream(hars = hars, genome_url_prefixes = GENOME_URL_PREFIXES).dump(output_file)
  output_file.close()

  shutil.copyfile('style.css', filejoin(output_dir, 'style.css'))

def print_summary(hars):
  for har in hars:
    for species in hars[har]['species']:
      print har, species
      img_dict = hars[har]['species'][species]['imgs']
      for img in sorted(img_dict.iterkeys()):
        print '  %2d: %s' % (img, img_dict[img])

def parse_args():
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('har_csv', help='Path to CSV file containing HAR IDs')
  parser.add_argument('img_dir', help='Path to top-level directory containing embryo images')
  parser.add_argument('--summary', action='store_true', help='Show list of images to be processed; will not produce a web page')
  parser.add_argument('--skipped', action='store_true', help='Show names of images that will not be processed')
  return parser.parse_args()

def main(args):
  csv_path = args.har_csv
  if not isfile(csv_path):
    print "Not a valid file: %s" % csv_path
    return
  hars = parse_har_csv(csv_path)
  img_dir_path = args.img_dir
  if not isdir(img_dir_path):
    print "Not a valid directory: %s" % img_dir_path
    return
  find_imgs(hars, img_dir_path, args.skipped)
  if args.summary:
    print_summary(hars)
    return

  output_dir = make_output_dir_name()
  os.mkdir(output_dir)
  setup_web_files(hars, output_dir)
  convert_imgs(hars, output_dir)

def ensure_pil():
  try:
    import Image
  except:
    print 'Could not import the Python Imaging library.'
    print 'To install it with pip: pip install PIL'
    sys.exit()

def ensure_jinja():
  try:
    import jinja2
  except:
    print 'Could not import the jinja2 library.'
    print 'To install it with pip: pip install jinja2'
    sys.exit()

def ensure_argparse():
  try:
    import argparse
  except:
    print 'Could not import the argparse library.'
    print 'This script requires Python 2.7 or greater'
    sys.exit()

ensure_argparse()
ensure_pil()
ensure_jinja()
main(parse_args())
