YUI.add("anim-node-plugin",function(a,b){var c=function(b){b=b?a.merge(b):{},b.node=b.host,c.superclass.constructor.apply(this,arguments)};c.NAME="nodefx",c.NS="fx",a.extend(c,a.Anim),a.namespace("Plugin"),a.Plugin.NodeFX=c},"3.14.1",{requires:["node-pluginhost","anim-base"]});
YUI.add("anim-scroll",function(a,b){var c=Number;a.Anim.behaviors.scroll={set:function(a,b,d,e,f,g,h){var i=a._node,j=[h(f,c(d[0]),c(e[0])-c(d[0]),g),h(f,c(d[1]),c(e[1])-c(d[1]),g)];j[0]&&i.set("scrollLeft",j[0]),j[1]&&i.set("scrollTop",j[1])},get:function(a){var b=a._node;return[b.get("scrollLeft"),b.get("scrollTop")]}}},"3.14.1",{requires:["anim-base"]});