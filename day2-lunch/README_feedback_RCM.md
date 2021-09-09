Looks good. Not a big deal, but it looks like you are getting extra spaces between lines in your output. Any ideas why and how you could avoid it?

Brendon: Fixed the extra lines. I forgot to strip the lines inputted from the mapper in remap.py, so I just had to add one .strip().
