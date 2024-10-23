import os
import yaml
import zipfile

with open('models/families.yml') as f:
	families = yaml.safe_load(f)

for (name, fam) in families.items():
	potts = fam['potts'].get('file', '')
	profile = fam['profile'].get('file', '')

	directory = os.path.dirname(potts)
	archive_file = os.path.join(directory, name + '_models.zip')
	with zipfile.ZipFile(archive_file, 'w', zipfile.ZIP_DEFLATED) as z:
		if os.path.isfile(potts):
			z.write(potts, os.path.basename(potts))
		else:
			print(f'File {potts} not found - skipping')

		if os.path.isfile(profile):
			z.write(profile, os.path.basename(profile))
		else:
			print(f'File {profile} not found - skipping')
