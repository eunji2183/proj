sudo apt install docker-io
docker --version
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu artful stable"
sudo docker run hello-world
sudo docker pull ensemblorg/ensembl-vep
sudo docker run -t -i ensemblorg/ensembl-vep ./vep
# Create a directory on your machine:
mkdir $HOME/vep_data

# Make sure that the created directory on your machine has read and write access granted
# so the docker container can write in the directory (VEP output):
chmod a+rwx $HOME/vep_data    

docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep
docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep perl INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 -g all
sudo docker run -t -i -v $HOME/vep_data:/opt/vep/.vep ensemblorg/ensembl-vep
cat /opt/vep/.vep/wes_cancer/project/config | while read id
do
echo "start vep_annotation for ${id} " `date`
./vep --cache --offline --format vcf --vcf --force_overwrite \
--dir_cache /opt/vep/.vep/ \
--dir_plugins /opt/vep/.vep/Plugins/ \
--input_file /opt/vep/.vep/wes_cancer/project/6.mutect/${id}_filter.vcf \
--output_file /opt/vep/.vep/wes_cancer/project/7.annotation/vep/${id}_vep.vcf 
echo "end vep_annotation for ${id} " `date`
done
