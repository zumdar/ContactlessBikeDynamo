# Contactless Bicycle Dynamo using Eddy Currents
![](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)
This is a small electrical generator that creates enough current to power a small bicycle light while your wheel is spinning. 


![alt text](./Documentation/dynamo_only.jpg)

![alt text](./Documentation/dynamo_and_light.PNG)


# Description
![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)


The design files here include the physical design files for the Dynamo and the LED light, as well as a custom PCB for the LED light which is a simple voltage rectifier doubler and choice of red or white led. 

If you would like to make this project:
- print off the files in:
  -  \Mechanical\Components\STLs\Final STL
- acquire parts in the ![Bill of Materials](./Documentation/Bill%20of%20Materials%20-%20Parts_REV1%20.csv)
- aquire tools in the ![Tool BoM](./Documentation/Bill%20of%20Materials%20-%20Tools%20Needed_REV1.csv)
- follow instructions in the ![ZINE](./Documentation/Workshop%20How%20To%20-%20bicycle%20dynamo%20zine.pdf)

This project started off with a desire to be a workshop I could host to show an interesting way to generate electricity and power a bike light without needing a battery. 

Here are the workshop slides:
https://docs.google.com/presentation/d/1txSlyZcSNcTIk4scu_-NPuIz-HN-7od_rSDEqWW2YtM/edit?usp=sharing


# Technical 
![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

## Physical Design Files
The parts are labelled as follows:

- P001-B100	-	Square Dynamo
- P001-B101	-	Dynamo Cover
- P001-B102	-	Magnet Holder Ring
- P001-B103	-	M5 Bearing Adaptor
- P001-B104	-	Fork Bracket
- P001-B105	-	PCB Cover
- P001-B106	-	LED Light Mount

There is also an assembly of the Dynamo labelled:

- P001-B001	-	Dynamo
## PCB Design files

The PCB contains red and white LEDs which can be selected by doing a solder bridge while assembling. 

You can't use both at the same time. If you solder both, only the Red LED will light (it has a lower forward voltage). 

The LED PCB design can be found in:
- \Electrical\Dyno_LED_PCB

![pcb render](./Documentation/PCB_3dView.jpg)


# Acknowledgments
![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)
Thank you to:
Lucas Ray, Clara Caspard, and YJ Rodrigues for helping develop this! 

# License
This is an open source project!

![](./Documentation/oshw_facts.svg)