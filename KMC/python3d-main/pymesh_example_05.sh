#sudo docker run -v $(pwd)/pymesh_examples:/pymesh_examples -it pymesh/pymesh bash -c 'pip install imageio;pip install matplotlib;python /pymesh_examples/pymesh_example_05.py'
sudo docker run -v $(pwd)/pymesh_examples:/pymesh_examples -it --rm pymesh/pymesh bash -c 'pip install imageio;pip install matplotlib;python /pymesh_examples/pymesh_example_05.py' | tee log.txt
