"use strict";YUI.add("rg-uri",function(a){a.namespace("rg");var b=function(){function a(a){return null!==a&&""!==a}function b(a){a=a||"";var b=/^(?:(?![^:@]+:[^:@\/]*@)([^:\/?#.]+):)?(?:\/\/)?((?:(([^:@]*)(?::([^:@]*))?)?@)?([^:\/?#]*)(?::(\d*))?)(((\/(?:[^?#](?![^?#\/]*\.[^?#\/.]+(?:[?#]|$)))*\/?)?([^?#\/]*))(?:\?([^#]*))?(?:#(.*))?)/,d=["source","protocol","authority","userInfo","user","password","host","port","relative","path","directory","file","query","anchor"],e={name:"queryKey",parser:/(?:^|&)([^&=]*)=?([^&]*)/g},f=b.exec(a),g=14,h=this;for(this.uriParts={};g--;)this.uriParts[d[g]]=f[g]||"";this.uriParts[e.name]={},this.uriParts[d[12]].replace(e.parser,function(a,b,c){b&&(h.uriParts[e.name][b]=c)}),this.queryObj=new c(this.uriParts.query),this.hasAuthorityPrefixUserPref=null}return b.prototype.protocol=function(a){return void 0!==a&&(this.uriParts.protocol=a),this.uriParts.protocol},b.prototype.hasAuthorityPrefix=function(a){return void 0!==a&&(this.hasAuthorityPrefixUserPref=a),null===this.hasAuthorityPrefixUserPref?-1!==this.uriParts.source.indexOf("//"):this.hasAuthorityPrefixUserPref},b.prototype.userInfo=function(a){return void 0!==a&&(this.uriParts.userInfo=a),this.uriParts.userInfo},b.prototype.host=function(a){return void 0!==a&&(this.uriParts.host=a),this.uriParts.host},b.prototype.port=function(a){return void 0!==a&&(this.uriParts.port=a),this.uriParts.port},b.prototype.path=function(a){return void 0!==a&&(this.uriParts.path=a),this.uriParts.path},b.prototype.query=function(a){return void 0!==a&&(this.queryObj=new c(a)),this.queryObj},b.prototype.anchor=function(a){return void 0!==a&&(this.uriParts.anchor=a),this.uriParts.anchor},b.prototype.setProtocol=function(a){return this.protocol(a),this},b.prototype.setHasAuthorityPrefix=function(a){return this.hasAuthorityPrefix(a),this},b.prototype.setUserInfo=function(a){return this.userInfo(a),this},b.prototype.setHost=function(a){return this.host(a),this},b.prototype.setPort=function(a){return this.port(a),this},b.prototype.setPath=function(a){return this.path(a),this},b.prototype.setQuery=function(a){return this.query(a),this},b.prototype.setAnchor=function(a){return this.anchor(a),this},b.prototype.getQueryParamValue=function(a){return this.query().getParamValue(a)},b.prototype.getQueryParamValues=function(a){return this.query().getParamValues(a)},b.prototype.deleteQueryParam=function(a,b){return 2===arguments.length?this.query().deleteParam(a,b):this.query().deleteParam(a),this},b.prototype.addQueryParam=function(a,b,c){return 3===arguments.length?this.query().addParam(a,b,c):this.query().addParam(a,b),this},b.prototype.replaceQueryParam=function(a,b,c){return 3===arguments.length?this.query().replaceParam(a,b,c):this.query().replaceParam(a,b),this},b.prototype.scheme=function(){var b="";return a(this.protocol())?(b+=this.protocol(),this.protocol().indexOf(":")!==this.protocol().length-1&&(b+=":"),b+="//"):this.hasAuthorityPrefix()&&a(this.host())&&(b+="//"),b},b.prototype.origin=function(){var b=this.scheme();return a(this.userInfo())&&a(this.host())&&(b+=this.userInfo(),this.userInfo().indexOf("@")!==this.userInfo().length-1&&(b+="@")),a(this.host())&&(b+=this.host(),a(this.port())&&(b+=":"+this.port())),b},b.prototype.toString=function(){var b=this.origin();return a(this.path())&&(b+=this.path()),a(this.query().toString())&&(0!==this.query().toString().indexOf("?")&&(b+="?"),b+=this.query().toString()),a(this.anchor())&&(0!==this.anchor().indexOf("#")&&(b+="#"),b+=this.anchor()),b},b.prototype.clone=function(){return new b(this.toString())},b}(),c=function(){function a(a){return a=decodeURIComponent(a),a=a.replace("+"," ")}function b(a){var b,c,d,e,f,g,h;if(this.params=[],void 0!==a&&null!==a&&""!==a)for(0===a.indexOf("?")&&(a=a.substring(1)),c=a.toString().split(/[&;]/),b=0;b<c.length;b++)d=c[b],e=d.split("="),g=e[0],f=null===e[1]?"":e[1],h=-1===d.indexOf("=")?null:f,this.params.push([g,h])}return b.prototype.getParamValue=function(b){var c,d;for(d=0;d<this.params.length;d++)if(c=this.params[d],a(b)===a(c[0]))return c[1]},b.prototype.getParamValues=function(b){var c,d,e=[];for(c=0;c<this.params.length;c++)d=this.params[c],a(b)===a(d[0])&&e.push(d[1]);return e},b.prototype.deleteParam=function(b,c){var d,e,f,g,h=[];for(d=0;d<this.params.length;d++)e=this.params[d],f=a(e[0])===a(b),g=a(e[1])===a(c),(1===arguments.length&&!f||2===arguments.length&&!f&&!g)&&h.push(e);return this.params=h,this},b.prototype.addParam=function(a,b,c){return 3===arguments.length&&-1!==c?(c=Math.min(c,this.params.length),this.params.splice(c,0,[a,b])):arguments.length>0&&this.params.push([a,b]),this},b.prototype.replaceParam=function(b,c,d){var e,f,g=-1;if(3===arguments.length){for(e=0;e<this.params.length;e++)if(f=this.params[e],a(f[0])===a(b)&&decodeURIComponent(f[1])===a(d)){g=e;break}this.deleteParam(b,d).addParam(b,c,g)}else{for(e=0;e<this.params.length;e++)if(f=this.params[e],a(f[0])===a(b)){g=e;break}this.deleteParam(b),this.addParam(b,c,g)}return this},b.prototype.toString=function(){var a,b,c="";for(a=0;a<this.params.length;a++)b=this.params[a],c.length>0&&(c+="&"),null===b[1]?c+=b[0]:c+=b.join("=");return c.length>0?"?"+c:c},b}();a.rg.Uri=b,a.rg.Query=c},"0.0.1",{});
"use strict";YUI.add("rg-ajax",function(a){a.namespace("rg");var b,c=0;a.rg.ajaxDbwAware=function(b,c,d,e,f,g,h){return a.Lang.isUndefined(g)&&(g="GET"),g=g.toUpperCase(),a.rg.ajax(b,c,d,e,f,g,h)},a.rg.ajaxEncodePostData=function(b,c){return function b(d,e,f){for(var g in e)if(e.hasOwnProperty(g)){var h=d?d+"["+g+"]":g;a.Lang.isObject(e[g])?b(h,e[g],f):void 0!==e[g]&&null!==e[g]&&(f[h]=c?e[g]:encodeURIComponent(e[g]))}return f}("",b,{})},a.rg.ajaxEncodePostDataToQueryString=function(b){var c=a.rg.ajaxEncodePostData(b),d=[];for(var e in c)c.hasOwnProperty(e)&&d.push(e+"="+c[e]);return d.length>0?d.join("&"):null};var d=function(d,e,f,g,h,i,j){var k=a.merge(j,{Accept:"application/json"});a.Lang.isObject(e)||(e={}),a.Lang.isUndefined(i)&&(i="POST"),"POST"!==i&&"PUT"!==i&&"DELETE"!==i&&"PATCH"!==i||(e.request_token=RGCommons.requestToken.get(),k["Content-Type"]="application/x-www-form-urlencoded"),b=i;var l;return{headers:k,method:i,data:a.rg.ajaxEncodePostDataToQueryString(e),arguments:{fn:f,context:g,args:h},on:{start:function(){l=(new Date).getTime()},success:function(b,k,m){var n,o=(new Date).getTime();if(204===k.status)n={success:!0,result:{},errors:[]};else try{if(RGCommons.performedRequests.log(d,i,k,!0,parseInt(o-l,10)),n=a.JSON.parse(k.responseText),n.debug&&a.Array.each(n.debug,function(b){a.log(b.message,"info"),a.log(b.context,"info")}),n.requestToken&&RGCommons.requestToken.update(n.requestToken),n.errors&&"csrf-retry"===n.errors[0]&&c<5)return c++,void a.rg.ajaxDbwAware(d,e,f,g,h,i,j);if(n.success&&(c=0),!1===n.success&&n.result.redirect)return a.log("ajax redirect: "+n.result.redirect),void a.rg.utils.url.forwardTo(n.result.redirect)}catch(a){n={success:!1,result:{},errors:["An error occurred while parsing server response.",a.message]}}a.Lang.isFunction(m.fn)&&m.fn.call(m.context,n,m.args)},failure:function(b,c,e){var f;try{RGCommons.performedRequests.log(d,i,c,!1),f=a.JSON.parse(c.responseText)}catch(a){f={success:!1,result:{},errors:["A server error occurred"]}}f.success=!1,!a.Lang.isFunction(e.fn)||c&&!c.status||e.fn.call(e.context,f,e.args)}}}};a.rg.ajax=function(c,e,f,g,h,i,j){if(!c)throw new Error("Trying to perform ajax request to empty url");if(b&&"GET"!==b||"GET"!==i){var k=new a.rg.Uri(c);k.addQueryParam("dbw","true"),c=k.toString()}return a.io(c,d(c,e,f,g,h,i,j))},a.rg.queueAjax=function(c,e,f,g,h,i,j){if(!c)throw new Error("Trying to queue ajax request to empty url");if(b&&"GET"!==b||"GET"!==i){var k=new a.rg.Uri(c);k.addQueryParam("dbw","true"),c=k.toString()}return a.io.queue(c,d(c,e,f,g,h,i,j))}},"0.0.1",{requires:["io-base","io-queue","json","cookie","rg-uri","rg-utils-url"]});
"use strict";YUI.add("rg-anim",function(a){a.namespace("rg"),a.rg.highlight=function(b,c,d){var e=a.one(b),f=a.merge({duration:1.5},d);if(0===a.UA.ie||a.UA.ie>8)e.setStyle("backgroundColor","#fffed2"),e.transition({easing:"ease-out",duration:f.duration,backgroundColor:"rgba(0,0,0,0)"},function(){e.setStyle("backgroundColor",""),c&&c()});else{var g=new a.Anim({node:e,from:{backgroundColor:"#fffed2"},to:{backgroundColor:"#ffffff"},duration:f.duration});g.set("easing",a.Easing.easeOut),g.on("end",function(){e.setStyle("backgroundColor",""),c&&c()}),g.run()}},a.rg.fadeToggle=function(b,c,d){var e,f,g=a.one(b);if(!g)return!1;e="none"!=g.getStyle("display"),f=a.merge({duration:.2},d),e?g.transition({easing:"linear",duration:f.duration,opacity:"0"},function(){g.hide(),g.setStyle("opacity",""),c&&c()}):(g.setStyle("opacity","0"),g.show(),g.transition({easing:"linear",duration:f.duration,opacity:"1"},function(){g.setStyle("opacity",""),c&&c()}))},a.rg.slideOpen=function(b,c,d){var e=a.one(b),f="none"!=e.getStyle("display"),g=null,h=d&&"duration"in d?d.duration:.3,i=d&&"dimension"in d?d.dimension:"height",j=new a.Anim({node:e,duration:h});if(!f){e.setStyle(i,""),e.setStyle("overflow","hidden"),e.setStyle("visibility","hidden"),e.setStyle("display","block");var k=e.getStyle("marginTop"),l=e.getStyle("marginBottom"),m=e.getStyle("marginLeft"),n=e.getStyle("marginRight"),o=e.getStyle("paddingTop"),p=e.getStyle("paddingBottom"),q=e.getStyle("paddingLeft"),r=e.getStyle("paddingRight");return"height"!=i?(e.setStyle("paddingLeft","0"),e.setStyle("paddingRight","0"),e.setStyle("marginLeft","0"),e.setStyle("marginRight","0"),g=e.get("offsetWidth"),j.set("to",{width:g,marginLeft:m,marginRight:n,paddingLeft:q,paddingRight:r})):(e.setStyle("paddingTop","0"),e.setStyle("paddingBottom","0"),e.setStyle("marginTop","0"),e.setStyle("marginBottom","0"),g=e.get("offsetHeight"),j.set("to",{height:g,marginTop:k,marginBottom:l,paddingTop:o,paddingBottom:p})),e.setStyle("display","none"),e.setStyle("visibility",""),e.setStyle(i,"1px"),j.on("start",function(){e.show(),e.setStyle("display","block")}),j.on("end",function(){e.removeAttribute("style"),c&&c(!0)}),j.run(),j}},a.rg.slideClose=function(b,c,d){var e=a.one(b);if(e){var f="none"!=e.getStyle("display"),g=d&&"duration"in d?d.duration:.2,h=d&&"dimension"in d?d.dimension:"height",i=new a.Anim({node:e,duration:g});if(f)return"height"!=h?i.set("to",{width:"1px",paddingLeft:"0",paddingRight:"0",marginLeft:"0",marginRight:"0"}):i.set("to",{height:"1px",paddingTop:"0",paddingBottom:"0",marginTop:"0",marginBottom:"0"}),i.on("start",function(){e.setStyle("overflow","hidden")}),i.on("end",function(){e.removeAttribute("style"),e.setStyle("display","none"),c&&c(!1)}),i.run(),i}},a.rg.slideToggle=function(b,c,d,e){var f=a.one(b);if(!f)return void(c&&c());var g="none"!=f.getStyle("display");return e||(e={}),e.dimension=d?"width":"height",g?a.rg.slideClose(b,c,e):a.rg.slideOpen(b,c,e)},a.rg.scrollTo=function(b,c,d){b=a.one(b);var e=new a.Anim({node:b,duration:d,to:{scroll:[0,c]},easing:a.Easing.easeInOut});e.run(),e.destroy()}},"0.0.1",{requires:["anim","transition"]});
"use strict";YUI.add("rg-base",function(){},"0.2.0",{requires:["array-extras","base","node","plugin","rg-ajax","event","querystring-stringify","lazysizes","rg-anim"]});