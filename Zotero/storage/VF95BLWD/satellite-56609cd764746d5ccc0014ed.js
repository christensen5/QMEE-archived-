var insertChunk = function(source, target) {

  var chunkData      = target + " data chunk",
      chunkTime      = target + " deposit time",
      request        = new XMLHttpRequest(),
      targetElement  = document.querySelector(target);

  request.open('GET', source, true);

  request.onload = function() {
    if (request.status >= 200 && request.status < 400 && targetElement) {

      var responseChunk = request.responseText,
          storedTime    = localStorage.getItem(chunkTime),
          storedChunk   = localStorage.getItem(chunkData),
          timeNow       = (new Date()).getTime(),
          oneDayAgo     = timeNow - (1000 * 60 * 60 * 24);

      if ( storedTime  === null || storedChunk === null || oneDayAgo > storedTime ) {

        localStorage.setItem(chunkData, responseChunk);
        localStorage.setItem(chunkTime, timeNow);

        targetElement.insertAdjacentHTML("afterbegin", responseChunk);
      } else {
        targetElement.insertAdjacentHTML("afterbegin", storedChunk);
      }

    } else {
      // We reached the server but couldn't find the source file or the target element is missing
      // Consider reporting errors via DTM
    }
  };

  request.onerror = function() {
    // There was a connection error while attempting to reach the server to fetch the source
  }
    request.send();
};
