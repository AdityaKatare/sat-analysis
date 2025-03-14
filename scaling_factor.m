uva_res = 240*1E-6 * 1.21;
irrad_uva = uva_res*(uva-6);

uvb_res = 144*1E-6 * 1.40;
irrad_uvb = uvb_res*(uvb-3);

uvc_res = 66.6*1E-6 * 0.90;
irrad_uvc_1hz = uvc_res*(uvc_1hz-5);
irrad_uvc = uvc_res*(uvc-5);