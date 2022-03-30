From debian
RUN echo "Hello from inside the container"
RUN apt-get -y update
RUN apt-get -y install git python3-tk python3-pip gfortran make texlive texlive-latex-extra texlive-fonts-recommended dvipng cm-super
RUN pip3 install numpy matplotlib

WORKDIR /app
COPY . .
RUN mkdir qgui
RUN mv * qgui; exit 0
RUN  pwd
WORKDIR /app/qgui
RUN pwd
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN python3 INSTALL.py root
RUN chmod a+x /usr/lib/QGUI/Qgui
WORKDIR /app
RUN rm -rf qgui



RUN git clone https://www.github.com/qusers/Q6.git 
WORKDIR /app/Q6/src
RUN make all
WORKDIR  /app/Q6/bin
RUN ls
RUN cp * /usr/bin
WORKDIR /app
RUN rm -rf Q6

CMD [ "/usr/bin/Qgui", "-p"]

