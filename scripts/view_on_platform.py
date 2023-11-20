def atoms_to_html(atoms):
    'Return the html representation the atoms object as string'
    from tempfile import NamedTemporaryFile
    with NamedTemporaryFile('r+', suffix='.html') as ntf:
        atoms.write(ntf.name, format='html')
        ntf.seek(0)
        html = ntf.read()
    return html

