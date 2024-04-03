sudo systemctl start docker
sudo docker build -t pore_designer:latest . --build-arg CACHEBUST=$(date +%s)