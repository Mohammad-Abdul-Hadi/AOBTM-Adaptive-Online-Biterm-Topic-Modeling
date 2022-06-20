FROM Ubuntu

RUN sudo apt update
RUN sudo apt upgrade
RUN sudo apt install python2.7
RUN sudo apt install python-pip

CMD ["sh" "./script/runExample.sh"]