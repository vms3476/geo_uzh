% script for testing coregistration of ALS and TLS data
if 0
  tls = load('/Users/morsdorf/Projects/Hyperforest/TLS/Wijn_30080.xyz');
  load /Users/morsdorf/Projects/Hyperforest/TLS/ALS_RAW_WIJN_30080.mat
else
  load /Users/morsdorf/Projects/Hyperforest/TLS/alstls_wijn_30080.mat
end
clf
myscatter3(raw.x,raw.y,raw.z,raw.z,ocean2(64),5);
hold on
myscatter3(tls(:,1),tls(:,2),tls(:,3),tls(:,3),gray(64),2);
mat3d2osg('alstls_wijn_30080.osg');