_satellite.pushBlockingScript(function(event, target, $variables){
  // Load SDI Data Layer Functions

if( typeof digitalData == 'undefined' ) { digitalData = {} }

_sdi.ddoSetVar = function(p, v) {
	for (var nodes = p.split("."), cursor = digitalData, i = 0; i < nodes.length - 1; i++)
		cursor[nodes[i]] || (cursor[nodes[i]] = {}),
            cursor = cursor[nodes[i]];
    return cursor[nodes[nodes.length-1]] = v,
	cursor[nodes[nodes.length-1]]
}

_sdi.ddoGetVar = function(p) {
  for (var cursor = digitalData, bits = p && "digitalData" !== p ? p.split(".") : "", i = 0; i < bits.length; i++) {
    if ("undefined" == typeof cursor[bits[i]]) {
      cursor = void 0;
      break
    }
    cursor = cursor[bits[i]]
  }
  return cursor;
}

_sdi.ddoDefault = function(p,v) {
  if( typeof _sdi.ddoGetVar(p) == "undefined" ) {
    _sdi.ddoSetVar(p,v);
  }
  return _sdi.ddoGetVar(p);
}

_sdi.ddoDefaults = function( _ddoDefaults ) {
  for( o in _ddoDefaults ) { 
    _sdi.ddoDefault(o,_ddoDefaults[o]);
  }
}
});
