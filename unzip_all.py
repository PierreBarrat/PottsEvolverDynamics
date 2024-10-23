import os
import yaml
import zipfile

with open('models/families.yml') as f:
	families = yaml.safe_load(f)

for (name, fam) in families.items():
	potts = fam['potts'].get('file', '')

	directory = os.path.dirname(potts)
	archive_file = os.path.join(directory, name + '_models.zip')
	if not os.path.isfile(archive_file):
		print(f'{archive_file} not found - skipping')
		
	with zipfile.ZipFile(archive_file, 'r') as archive:
		archive.extractall(directory)
