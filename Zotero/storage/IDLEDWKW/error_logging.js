!function(a){function b(a,b){try{if("function"!=typeof a)return a;if(!a.bugsnag){var c=d();a.bugsnag=function(d){if(b&&b.eventHandler&&(r=d),s=c,!v){var e=a.apply(this,arguments);return s=null,e}try{return a.apply(this,arguments)}catch(a){throw u.notifyException(a),p(),a}finally{s=null}},a.bugsnag.bugsnag=a.bugsnag}return a.bugsnag}catch(b){return a}}function c(){y=!1}function d(){var a=document.currentScript||s;if(!a&&y){var b=document.scripts||document.getElementsByTagName("script");a=b[b.length-1]}return a}function e(a){var b=d();b&&(a.script={src:b.src,content:b.innerHTML?b.innerHTML.substr(0,50):void 0})}function f(b){var c=a.console;void 0!==c&&void 0!==c.log&&c.log("[Bugsnag] "+b)}function g(b,c,d){if(d>=x)return encodeURIComponent(c)+"=[RECURSIVE]";d=d+1||1;try{if(a.Node&&b instanceof a.Node)return encodeURIComponent(c)+"="+encodeURIComponent(o(b));var e=[];for(var f in b)if(b.hasOwnProperty(f)&&null!=f&&null!=b[f]){var h=c?c+"["+f+"]":f,i=b[f];e.push("object"==typeof i?g(i,h,d):encodeURIComponent(h)+"="+encodeURIComponent(i))}return e.sort().join("&")}catch(a){return encodeURIComponent(c)+"="+encodeURIComponent(""+a)}}function h(a,b){var c=new XMLHttpRequest;c.open("post",a),c.setRequestHeader("Content-Type","application/json"),c.setRequestHeader("Accept","application/json"),c.send(g(b))}function i(b,c){var d=a.rgConfig&&void 0!==a.rgConfig[b]?a.rgConfig[b]:c;return"false"===d&&(d=!1),d}function j(b,c){if(a.atob&&!k()&&("undefined"!=typeof RGCommons&&RGCommons.debug&&RGCommons.debug.logException({Type:"JavascriptException",PreviousType:b.name,message:b.message,File:b.file,Line:b.lineNumber+":"+b.columnNumber}),!b.file||!b.file.match(/\/pdf\.viewer\.js$/))){var d=[b.name,b.message,b.stacktrace].join("|");if(d!==t){t=d,c=c||{},"object"!=typeof c&&(c={metaData:c}),r&&(c["Last Event"]=n(r));var e={AccountId:i("accountId"),Module:i("module"),Action:i("action"),PageId:i("pageId"),CorrelationId:i("correlationId"),jsMetaData:c,Request:a.location.href,UserAgent:navigator.userAgent,Type:"JavascriptException",PreviousType:b.name,message:b.message,Stacktrace:b.stacktrace,File:b.file,Line:b.lineNumber,Column:b.columnNumber};if(0===e.Line&&/Script error\.?/.test(e.message))return f("Ignoring cross-domain script error.");h("go.Error.html",e)}}}function k(){return!!a.frameElement}function l(){var a;try{throw new Error("")}catch(b){a=m(b)}if(!a){var b=[];try{for(var c=arguments.callee.caller.caller;c&&b.length<10;){var d=z.test(c.toString())?RegExp.$1||"[anonymous]":"[anonymous]";b.push(d),c=c.caller}}catch(a){f(a)}a=b.join("\n")}return a}function m(a){return a.stack||a.backtrace||a.stacktrace}function n(a){return{millisecondsAgo:new Date-a.timeStamp,type:a.type,which:a.which,target:o(a.target)}}function o(a){if(a){var b=a.attributes;if(b){for(var c="<"+a.nodeName.toLowerCase(),d=0;d<b.length;d++)b[d].value&&"null"!==b[d].value.toString()&&(c+=" "+b[d].name+'="'+b[d].value+'"');return c+">"}return a.nodeName}}function p(){w+=1,a.setTimeout(function(){w-=1})}function q(a,b,c){var d=a[b];a[b]=c(d)}var r,s,t,u={},v=!0,w=0,x=5;u.notifyException=function(a,b,c){a&&(b&&"string"!=typeof b&&(c=b,b=void 0),c||(c={}),e(c),j({name:b||a.name,message:a.message||a.description,stacktrace:m(a)||l(),file:a.fileName||a.sourceURL,lineNumber:a.lineNumber||a.line,columnNumber:a.columnNumber?a.columnNumber+1:void 0},c))},u.notify=function(b,c,d){b||(b="RGNotify",c="Bugsnag.notify() was called with no arguments",f(c)),j({name:b,message:c,stacktrace:l(),file:a.location.toString(),lineNumber:1},d)};var y="complete"!==document.readyState;document.addEventListener?(document.addEventListener("DOMContentLoaded",c,!0),a.addEventListener("load",c,!0)):a.attachEvent("onload",c);var z=/function\s*([\w\-$]+)?\s*\(/i;if(a.atob){if(a.ErrorEvent)try{0===new a.ErrorEvent("test").colno&&(v=!1)}catch(a){}}else v=!1;q(a,"onerror",function(b){return function(c,d,f,g,h){var i={};!g&&a.event&&(g=a.event.errorCharacter),e(i),s=null,w||("string"!=typeof c&&(c=JSON.stringify(c)),j({name:h&&h.name||"window.onerror",message:c,file:d,lineNumber:f,columnNumber:g,stacktrace:h&&m(h)||l()},i)),b&&b(c,d,f,g,h)}});var A=function(a){return function(c,d){if("function"==typeof c){c=b(c);var e=Array.prototype.slice.call(arguments,2);return a(function(){c.apply(this,e)},d)}return a(c,d)}};q(a,"setTimeout",A),q(a,"setInterval",A),a.requestAnimationFrame&&q(a,"requestAnimationFrame",function(a){return function(c){return a(b(c))}}),a.setImmediate&&q(a,"setImmediate",function(a){return function(){var c=Array.prototype.slice.call(arguments);return c[0]=b(c[0]),a.apply(this,c)}}),"onunhandledrejection"in a&&a.addEventListener("unhandledrejection",function(a){var b=a.reason;b&&(b instanceof Error||b.message)?u.notifyException(b):u.notify("UnhandledRejection",b)}),"EventTarget Window Node ApplicationCache AudioTrackList ChannelMergerNode CryptoOperation EventSource FileReader HTMLUnknownElement IDBDatabase IDBRequest IDBTransaction KeyOperation MediaController MessagePort ModalWindow Notification SVGElementInstance Screen TextTrack TextTrackCue TextTrackList WebSocket WebSocketWorker Worker XMLHttpRequest XMLHttpRequestEventTarget XMLHttpRequestUpload".replace(/\w+/g,function(c){var d=a[c]&&a[c].prototype;d&&d.hasOwnProperty&&d.hasOwnProperty("addEventListener")&&(q(d,"addEventListener",function(a){return function(c,d,e,g){try{d&&d.handleEvent&&(d.handleEvent=b(d.handleEvent,{eventHandler:!0}))}catch(a){f(a)}return a.call(this,c,b(d,{eventHandler:!0}),e,g)}}),q(d,"removeEventListener",function(a){return function(c,d,e,f){return a.call(this,c,d,e,f),a.call(this,c,b(d),e,f)}}))}),a.Bugsnag=u}(window);