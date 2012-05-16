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
label.rt(btex $|P|=L/10$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw dot shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=1;
label.rt(btex $|P|=L$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
%draw cross shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

ls:=2;
label.rt(btex $|P|=5$ etex, (1cm,sh + ls*ups));
draw (sx,sh + ls*ups)--(1cm,sh + ls*ups);
draw circle shifted 1/2[(sx,sh + ls*ups),(1cm,sh + ls*ups)];

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

  gdraw "MackenthunCoherentM4L1000absP100a"; gdraw "MackenthunCoherentM4L1000absP100a" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "MackenthunCoherentM4L100absP10a"; gdraw "MackenthunCoherentM4L100absP10a" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "MackenthunCoherentM4L10absP1a"; gdraw "MackenthunCoherentM4L10absP1a" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunCoherentM4L10000absP1000a"; gdraw "MackenthunCoherentM4L10000absP1000a" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "MackenthunCoherentM4L1000absP5a"; gdraw "MackenthunCoherentM4L1000absP5a" plot(circle); 
  gdraw "MackenthunCoherentM4L100absP5a"; gdraw "MackenthunCoherentM4L100absP5a" plot(circle); 
  gdraw "MackenthunCoherentM4L10absP5a"; gdraw "MackenthunCoherentM4L10absP5a" plot(circle);

  gdraw "MackenthunCoherentM4L1000absP1000a"; 
  gdraw "MackenthunCoherentM4L100absP100a"; 
  gdraw "MackenthunCoherentM4L10absP10a"; 


draw legend shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of phase error

beginfig(2);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE phase etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,"1.1e-6", whatever,90);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "cltM4L1000absP100p"; gdraw "MackenthunCoherentM4L1000absP100p" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "cltM4L100absP10p"; gdraw "MackenthunCoherentM4L100absP10p" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "cltM4L10absP1p"; gdraw "MackenthunCoherentM4L10absP1p" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunCoherentM4L10000absP1000p"; gdraw "MackenthunCoherentM4L10000absP1000p" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "cltM4L1000absP5p"; gdraw "MackenthunCoherentM4L1000absP5p" plot(circle); 
  gdraw "cltM4L100absP5p"; gdraw "MackenthunCoherentM4L100absP5p" plot(circle); 
  gdraw "cltM4L10absP5p"; gdraw "MackenthunCoherentM4L10absP5p" plot(circle);

  gdraw "MackenthunCoherentM4L1000absP1000p"; 
  gdraw "MackenthunCoherentM4L100absP100p"; 
  gdraw "MackenthunCoherentM4L10absP10p";

draw legend shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;


%%%%%%%%%%%%
%%%% Begin plot of complex error

beginfig(3);
draw begingraph(plotwidth, plotheight);
 setcoords(linear,log);
 glabel.bot(btex SNR (dB) etex, OUT) shifted (0cm,1.3cm);
  glabel.lft(btex MSE complex gain etex, OUT) rotated 90 shifted (1.5cm,1cm);
  
  %setrange(whatever,"1.1e-6", whatever,90);
  %setrange(4.9,"1e-13", 22,2);

%  gdraw (0, "1e-13")--(0, 2) dashed evenly; %glabel.bot(btex AWGN channel capacity etex, 0) rotated 90;

  gdraw "MackenthunCoherentM4L1000absP100c"; gdraw "MackenthunCoherentM4L1000absP100c" plot(dot); glabel.llft(btex $L = 1000$ etex, 1000);
  gdraw "MackenthunCoherentM4L100absP10c"; gdraw "MackenthunCoherentM4L100absP10c" plot(dot); glabel.llft(btex $L = 100$ etex, 1000);
  gdraw "MackenthunCoherentM4L10absP1c"; gdraw "MackenthunCoherentM4L10absP1c" plot(dot); glabel.llft(btex $L = 10$ etex, 1000);
  %gdraw "MackenthunCoherentM4L10000absP1000c"; gdraw "MackenthunCoherentM4L10000absP1000c" plot(dot); glabel.llft(btex $L = 10000$ etex, 1000);

  gdraw "MackenthunCoherentM4L1000absP5c"; gdraw "MackenthunCoherentM4L1000absP5c" plot(circle);
  gdraw "MackenthunCoherentM4L100absP5c"; gdraw "MackenthunCoherentM4L100absP5c" plot(circle);
  gdraw "MackenthunCoherentM4L10absP5c"; gdraw "MackenthunCoherentM4L10absP5c" plot(circle);

  gdraw "MackenthunCoherentM4L1000absP1000c"; 
  gdraw "MackenthunCoherentM4L100absP100c"; 
  gdraw "MackenthunCoherentM4L10absP10c";

draw legend shifted(legendx,legendy);

endgraph;

%clip currentpicture to (-10cm,-10cm)--(-10cm,plotheight)--(plotwidth+5cm,plotheight)--(plotwidth+5cm,-10cm)--cycle

endfig;

end;