# Design Symmetric Pores
Use a workflow including ProteinMPNN and AlphaFold2 to optimize the sequence of an barrel/pore from a starting scaffold.

# Running
Place the scaffold PDB in a folder and then:

1. Ensure docker deamon is started `sudo systemctl start docker`

2. Interactively start the image `sudo docker run -it --entrypoint /bin/bash -v $(pwd)/:/workspace --gpus all antiquatedarachnid/pore_designer:latest`

3. Setup a config file for the design run `python ../pore_designer/scripts/design_pore.py make-config output input.pdb`
- use `--help` to see additional options that can be given to alter stringency of the run and give better sequences.

4. Start the run `python ../pore_designer/scripts/design_pore.py design-pore output/config.json`
