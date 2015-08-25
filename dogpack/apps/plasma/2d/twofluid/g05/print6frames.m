function print6frames(n)

sn = num2str(n);

% File names
fname01 = 'rhoi_xx_color.jpg';
if n<10
 fname01(6)='0';
 fname01(7)=sn(1);
else
 fname01(6)=sn(1);
 fname01(7)=sn(2);
end

fname03 = 'rhoe_xx_color.jpg';
if n<10
 fname03(6)='0';
 fname03(7)=sn(1);
else
 fname03(6)=sn(1);
 fname03(7)=sn(2);
end

fname05 = 'veli_xx_color.jpg';
if n<10
 fname05(6)='0';
 fname05(7)=sn(1);
else
 fname05(6)=sn(1);
 fname05(7)=sn(2);
end

fname07 = 'vele_xx_color.jpg';
if n<10
 fname07(6)='0';
 fname07(7)=sn(1);
else
 fname07(6)=sn(1);
 fname07(7)=sn(2);
end

fname09 = 'Bmag_xx_color.jpg';
if n<10
 fname09(6)='0';
 fname09(7)=sn(1);
else
 fname09(6)=sn(1);
 fname09(7)=sn(2);
end

fname11 = 'Emag_xx_color.jpg';
if n<10
 fname11(6)='0';
 fname11(7)=sn(1);
else
 fname11(6)=sn(1);
 fname11(7)=sn(2);
end

disp([' fname01 = ',fname01]);
disp([' fname03 = ',fname03]);
disp([' fname05 = ',fname05]);
disp([' fname07 = ',fname07]);
disp([' fname09 = ',fname09]);
disp([' fname11 = ',fname11]);

figure(1);  print('-djpeg95',fname01);
figure(3);  print('-djpeg95',fname03);
figure(5);  print('-djpeg95',fname05);
figure(7);  print('-djpeg95',fname07);
figure(9);  print('-djpeg95',fname09);
figure(11); print('-djpeg95',fname11);
