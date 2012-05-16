prologues := 2;
filenametemplate "%j-%c.mps";

input TEX; 
TEXPRE("%&latex" & char(10) & 
"\documentclass{minimal}\begin{document}"); 
TEXPOST("\end{document}"); 

input graph;
input sarith;
picture ORIGINAL_PICTURE;
ORIGINAL_PICTURE := currentpicture;

picture dd; dd :=
dashpattern(on 4bp off 2bp on 1bp off 2bp);

picture diamond;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm rotated 45;
diamond := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture square;
draw unitsquare shifted(-0.5,-0.5) scaled 0.17cm;
square := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture circle; 
draw fullcircle scaled 0.19cm;
circle := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture cross; 
draw (-0.09cm,0.0)--(0.09cm,0.0);
draw (0,-0.09cm)--(0,0.09cm);
cross := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture dot; 
fill fullcircle scaled 0.15cm;
dot := currentpicture; currentpicture := ORIGINAL_PICTURE;

picture times; times := cross rotated 45;


%%%%%%%%%%%%%
%%%%%Legend picture

picture legend;
sx = 0.2cm;
sh = 0.2cm;
ups = 0.35cm;

ls:=0;
label.rt(btex $|P|=L/4$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex $|P|=L/8$ etex, (1cm,sh + ls*ups));
%draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=3;
label.rt(btex Theory etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];


legend := currentpicture; currentpicture := ORIGINAL_PICTURE;


plotheight = 13cm;
plotwidth = 10cm;

legendx = 0.1cm;
legendy = 0.1cm;


%%%%%%%%%%%%
%%%% Begin plot of amplitude error

beginfig(1);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE amplitude etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,whatever, whatever,100);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "MackenthunCoherentM2L2048absP512a"; gdraw "MackenthunCoherentM2L2048absP512a" plot(dot); glabel.llft(btex $L = 2048$ etex, 1000);
  gdraw "MackenthunCoherentM2L256absP64a"; gdraw "MackenthunCoherentM2L256absP64a" plot(dot); glabel.llft(btex $L = 256$ etex, 1000);
  gdraw "MackenthunCoherentM2L32absP8a"; gdraw "MackenthunCoherentM2L32absP8a" plot(dot); glabel.llft(btex $L = 32$ etex, 1000);
  %gdraw "MackenthunCoherentM2L10000absP1000a"; gdraw "MackenthunCoherentM2L10000absP1000a" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "MackenthunCoherentM2L2048absP128a"; gdraw "MackenthunCoherentM2L2048absP128a" plot(circle); 
  gdraw "MackenthunCoherentM2L256absP16a"; gdraw "MackenthunCoherentM2L256absP16a" plot(circle); 
  gdraw "MackenthunCoherentM2L32absP2a"; gdraw "MackenthunCoherentM2L32absP2a" plot(circle);

  gdraw "MackenthunCoherentM2L2048absP2048a"; 
  gdraw "MackenthunCoherentM2L256absP256a"; 
  gdraw "MackenthunCoherentM2L32absP32a"; 


draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of phase error

beginfig(2);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE phase etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  setrange(whatever,whatever, whatever,5);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "cltM2L2048absP512p"; gdraw "MackenthunCoherentM2L2048absP512p" plot(dot); glabel.llft(btex $L = 2048$ etex, 1000);
  gdraw "cltM2L256absP64p"; gdraw "MackenthunCoherentM2L256absP64p" plot(dot); glabel.llft(btex $L = 256$ etex, 1000);
  gdraw "cltM2L32absP8p"; gdraw "MackenthunCoherentM2L32absP8p" plot(dot); glabel.llft(btex $L = 32$ etex, 1000);
  %gdraw "MackenthunCoherentM2L2560absP1000p"; gdraw "MackenthunCoherentM2L2560absP1000p" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "cltM2L2048absP128p"; gdraw "MackenthunCoherentM2L2048absP128p" plot(circle); 
  gdraw "cltM2L256absP16p"; gdraw "MackenthunCoherentM2L256absP16p" plot(circle); 
  gdraw "cltM2L32absP2p"; gdraw "MackenthunCoherentM2L32absP2p" plot(circle);

  gdraw "cltM2L2048absP2048p"; gdraw "MackenthunCoherentM2L2048absP2048p" plot(cross); 
  gdraw "cltM2L256absP256p"; gdraw "MackenthunCoherentM2L256absP256p" plot(cross); 
  gdraw "cltM2L32absP32p"; gdraw "MackenthunCoherentM2L32absP32p" plot(cross); 

draw legend shifted(legendx,legendy);

endgraph;

clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;

end;