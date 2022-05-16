y1 = board_adc_data;
d1 = amplifier_data;

DD2 = DD2;
Y2 = Y2;


y1 = board_adc_data;
d1 = amplifier_data;

y2 = board_adc_data;
d2 = amplifier_data;

y3 = board_adc_data;
d3 = amplifier_data;

y4 = board_adc_data;
d4 = amplifier_data;




DD2 = [d1 d2 d3 d4];
Y2 = [y1 y2 y3 y4];

DD2 = [d1 d2 d3];
Y2 = [y1 y2 y3];


y1 = board_adc_data;
d1 = amplifier_data;

y2 = board_adc_data;
d2 = amplifier_data;


DD2 = [d1 d2];
Y2 = [y1 y2];


save('DD2.mat','DD2')
save('Y2.mat','Y2')
