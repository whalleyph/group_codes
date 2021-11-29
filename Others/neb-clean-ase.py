echo -n "This will clean the current ASE-NEB calculation. Are you sure? (y: yes): "
read ans

if [ "$ans" != "y" ]; then
echo "Nothing is done..."
exit 1
fi

rm -rf IMG* image-* idpp.* neb.* slurm-* machine.file.* diffusion-barrier.png
