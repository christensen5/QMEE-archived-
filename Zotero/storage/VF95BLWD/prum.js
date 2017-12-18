!function(){"use strict";function e(e,t){return Object.keys(e).map(function(n){return n+"="+(t?e[n]:encodeURIComponent(e[n]))}).join("&")}function t(e,t){var n,o,r={};return e&&"string"==typeof e?(e.trim().split("&").forEach(function(e){n=e.indexOf("="),o=e.substring(n+1),r[e.substr(0,n)]=t?o:decodeURIComponent(o)}),r):r}function n(e){return e instanceof Date&&(e=e.valueOf()),"number"==typeof e&&parseInt(e.toString().substring(0,10),10)}function o(){var e={},t=arguments;t[0]instanceof Array&&(t=t[0]);for(var n=0,s=t.length;n<s;n++){var a=t[n];if(r(a))a=o(a);else if(!i(a))continue;Object.keys(a).forEach(function(t){a.hasOwnProperty(t)&&(e[t]=a[t])})}return e}function r(e){return"[object Array]"===Object.prototype.toString.call(e)}function i(e){return"[object Object]"===Object.prototype.toString.call(e)}function s(t,n,r,s){function a(t){return r+"?"+e(t)}function c(c,u){var d=s||!1;if(i(u)&&(u=o({id:n.getSiteID()},n.getSessionInfo(),u)),"GET"===c&&(r=a(u)),-1!==t.navigator.appName.indexOf("Internet Explorer")){var f=t.navigator.appVersion.match(/MSIE (\d+)/);f&&parseInt(f[1])<=9&&(d=!0)}if(t.XMLHttpRequest&&!d){var g=new t.XMLHttpRequest;g.open(c,r),"GET"===c?g.send(e(u)):g.send(JSON.stringify(u))}else t.document.createElement("img").src=r}return{get:function(e){c("GET",e)},post:function(e){c("POST",e)}}}function a(n,o){function r(e){return o.storageKey+"="+e+"; expires="+i}var i=new Date(Date.now()+o.retVisitor),s={usesCookies:!0,getItem:function(e){return s._cookieToObject()[e]},setItem:function(e,t){var n=s._cookieToObject(),o=n||{};o[e]=t,s._objectToCookie(o)},removeItem:function(e){var t=s._cookieToObject();t&&(delete t[e],s._objectToCookie(t))},_cookieToObject:function(){return t(n.document.cookie.split(";").filter(function(e){return e.length&&e.indexOf(o.storageKey)>-1}).join("").replace(o.storageKey+"=",""))},_objectToCookie:function(t){n.document.cookie=r(e(t))}};return s}function c(n,o){function r(e){if(c&&e){var n=t(c.getItem(o.storageKey));return n?n[e]:""}return""}function i(n,r){if(c&&n)try{var i=c.getItem(o.storageKey),s=i?t(i):{};s[n]=r,c.setItem(o.storageKey,e(s))}catch(e){console.error("unable to store "+n+" in storage.",e)}}function s(e){c&&e&&c.removeItem(e)}var c;return function(){var e=o.storageKey+"_enabled";if(n.localStorage&&"1"===n.localStorage.getItem(e))return void(c=n.localStorage);if(n.localStorage)try{if(n.localStorage.setItem(e,1),"1"===n.localStorage.getItem(e))return void(c=n.localStorage)}catch(e){console.error("localStorage.setItem() failed. Using cookies")}c=a(n,o)}(),{get:r,set:i,remove:s}}function u(e,t,o){var i={modules:[],send:void 0,storage:void 0,storageKey:e.storageKey,id:e.id,url:e.url,rum1url:e.rum1url,ver:e.ver,sessionIDLength:e.sessionIDLength,sessionLifetime:e.sessionLifetime,retVisitor:e.retVisitor,getSessionInfo:function(){var e=i.storage.get("sid"),t=parseInt(i.storage.get("sst"),10),o=n(Date.now());return e&&t?o-t>i.sessionLifetime?i.sessionStart(o-t<i.retVisitor):{sId:e,sST:t,sIS:i.getSessionInteractionStep(),rV:i.storage.get("rv")||"0",v:i.ver}:i.sessionStart(!1)},generateSessionID:function(){return(78364164096+Math.floor(2742745743359*Math.random())).toString(36)},sessionStart:function(e){e=e?"1":"0";var t=i.generateSessionID();i.storage.set("sid",t);var n=i.sessionMarkActive();return i.storage.set("sis","1"),i.storage.set("rv",e),{sId:t,sST:n,sIS:"1",rV:e,v:i.ver}},sessionMarkActive:function(){var e=n(Date.now());return i.storage.set("sst",e),e},getSessionInteractionStep:function(){return parseInt(i.storage.get("sis"),10)||1},bumpSessionInteractionStep:function(){i.storage.set("sis",i.getSessionInteractionStep()+1)},checkBrowser:function(){return t.document&&t.document.readyState&&Array.prototype.forEach&&Array.prototype.map},getSiteID:function(){return i.id.length||r(t._prum)&&r(t._prum[0])&&"id"===t._prum[0][0]&&(i.id=t._prum[0][1],i.storage.set("r1","1")),i.id},initialize:function(){i.storage=c(t,i),i.getSessionInfo(),i.send=s(t,i,i.url),i.sendLegacy=s(t,i,i.rum1url,!0),r(o)&&(i.modules=o,i.modules.forEach(function(e){e(t,i)}))}};i.checkBrowser()&&("complete"!==t.document.readyState?t.addEventListener("load",function e(t){t.target.removeEventListener("load",e),i.initialize()}):i.initialize())}function d(e){return{sAW:e.screen.availWidth,sAH:e.screen.availHeight,bIW:e.innerWidth,bIH:e.innerHeight,pD:e.screen.pixelDepth,dPR:1|e.devicePixelRatio,or:e.screen.orientation&&e.screen.orientation.type||""}}function f(e,t){function n(e){return e>0?e-u.navigationStart:-1}function r(n){var r=e.location;n.push({s:"nt",title:e.document.title,path:r.protocol+"//"+r.host+r.pathname,ref:e.document.referrer});var i=o(n);1!==i.tpe&&(t.send.get(i),"1"===t.storage.get("r1")&&t.sendLegacy.get(i),t.bumpSessionInteractionStep())}function i(){return{nT:f.navigation.type,rC:f.navigation.redirectCount}}function s(){return"https:"===e.location.protocol&&u.secureConnectionStart>0?n(u.secureConnectionStart):-1}function a(){return{nS:0,cS:n(u.connectStart),cE:n(u.connectEnd),dLE:n(u.domainLookupEnd),dLS:n(u.domainLookupStart),fS:n(u.fetchStart),hS:s(),rE:n(u.redirectEnd),rS:n(u.redirectStart),reS:n(u.requestStart),resS:n(u.responseStart),resE:n(u.responseEnd),uEE:n(u.unloadEventEnd),uES:n(u.unloadEventStart),dL:n(u.domLoading),dI:n(u.domInteractive),dCLES:n(u.domContentLoadedEventStart),dCLEE:n(u.domContentLoadedEventEnd),dC:n(u.domComplete),lES:n(u.loadEventStart),lEE:n(u.loadEventEnd)}}var c,u,f=e.performance||{};!function(){(u=f.timing)&&(c=setTimeout(function(){if(f.timing.loadEventEnd){clearInterval(c);var t=[];t.push(d(e)),t.push(i()),t.push(a()),r(t)}},25))}()}!function(e){var t=[f];u({storageKey:"pa-l",id:"",url:"//rum-collector-2.pingdom.net/img/beacon.gif",rum1url:"//rum-collector.pingdom.net/img/beacon.gif",ver:"1.3.0",sessionIDLength:parseInt("8",10),sessionLifetime:parseInt("1800",10),retVisitor:24*parseInt("30",10)*3600},e,t)}(window)}();