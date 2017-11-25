(function ($) {

  /**
   * Helper function to check browser support
   * @return
   *   Object containing
   *   -- browser:        Name of the browser
   *   -- browserVersion: Version of the browser
   *   -- acrobat:        Bit to check acrobat support
   *   -- acrobatVersion: Version of acrobat if present
   */
  var getAcrobatInfo = function() {

    var getBrowserName = function() {
        var userAgent = navigator ? navigator.userAgent.toLowerCase() : "other";
        if(userAgent.indexOf("chrome") > -1)        return "chrome";
        else if(userAgent.indexOf("safari") > -1)   return "safari";
        else if(userAgent.indexOf("msie") > -1)     return "ie";
        else if(userAgent.indexOf("firefox") > -1)  return "firefox";
        else if(userAgent.match(/trident.*rv[ :]*11\./))  return "ie"; // this is fix for IE 11
        return userAgent;
    };

    var getBrowserVersion = function() {
        var ua = navigator.userAgent, temp;
        // this is fix for IE 11
        if (ua.toLowerCase().match(/trident.*rv[ :]*11\./)) {
          return 11;
        }
        var M = ua.match(/(opera|chrome|safari|firefox|msie)\/?\s*(\.?\d+(\.\d+)*)/i);
        if(M && (temp= ua.match(/version\/([\.\d]+)/i))!= null) {
          M[2]= temp[1];
        }
        browserVersion = M ? M[2] : navigator.appVersion;
        return browserVersion;
    };

    var getActiveXObject = function(name) {
      try { return new ActiveXObject(name); } catch(e) {}
    };

    var getNavigatorPlugin = function(name) {
      for(key in navigator.plugins) {
        var plugin = navigator.plugins[key];
        if(plugin.name == name) return plugin;
      }
    };

    var getPDFPlugin = function() {
      return this.plugin = this.plugin || function() {
        if(getBrowserName() == 'ie') {
          //
          // load the activeX control
          // AcroPDF.PDF is used by version 7 and later
          // PDF.PdfCtrl is used by version 6 and earlier
          return getActiveXObject('AcroPDF.PDF') || getActiveXObject('PDF.PdfCtrl');
        }
        else {
          return getNavigatorPlugin('Adobe Acrobat') || getNavigatorPlugin('Chrome PDF Viewer') || getNavigatorPlugin('WebKit built-in PDF');
        }
      }();
    };

    var isAcrobatInstalled = function() {
      return !!getPDFPlugin();
    };

    var getAcrobatVersion = function() {
      try {
        var plugin = getPDFPlugin();

        if(getBrowserName() == 'ie') {
          var versions = plugin.GetVersions().split(',');
          var latest   = versions[0].split('=');
          return parseFloat(latest[1]);
        }

        if(plugin.version) return parseInt(plugin.version);
        return plugin.name

      }
      catch(e) {
        return null;
      }
    }

    //
    // The returned object
    //
    return {
      browser:        getBrowserName(),
      browserVersion: getBrowserVersion(),
      acrobat:        isAcrobatInstalled() ? 'installed' : false,
      acrobatVersion: getAcrobatVersion()
    };
  };

  Drupal.behaviors.highwire_toggle_pdf = {
    attach: function (context, settings) {
      // get the brower parameters
      var browserInfo = getAcrobatInfo();
      var browser_name = browserInfo.browser;
      var browser_version = browserInfo.browserVersion;
      var acrobat = browserInfo.acrobat;
      var acrobatVersion = browserInfo.acrobatVersion;

      var pdf_path = $('div.highwire-pdf').attr('data-pdf-path');
      var pdfjs_path = $('div.highwire-pdf').data('pdfjs-path');
      var pdf_container = '';

      if ((browser_name == "chrome") || (browser_name == "firefox") || ((browser_name == "safari") && (parseInt(browser_version) >= 5.0))) {
        show_pdf = true;
        if (acrobat) {
          var pdf_container = '<object width="100%"><embed src="' + pdf_path + '#zoom=75" width="100%" height="785px"></object>';
        }
        else {
          pdf_container =  "<iframe src='" + pdfjs_path + "' width='100%' height='750px'></iframe>";
        }

        $('.highwire-pdf').html(pdf_container);
        $('.higwire-pdf-download-link').remove();
      }
      else if(((browser_name == "safari") && (parseInt(browser_version) < 5.0))) {
        show_pdf = true;
        if (acrobat) {
          $('.higwire-pdf-download-link').remove();

          var pdf_container = '<object width="100%"><embed src="' + pdf_path + '#zoom=80" width="100%" height="785px"></object>';
        }
        else{
          pdf_container = '<iframe src="http://docs.google.com/viewer?url=' + pdf_path +'&embedded=true" style="width:718px; height:700px;" frameborder="0"></iframe>';
        }
        $('.highwire-pdf').html(pdf_container);
      }
      else if ((browser_name == "ie") && (parseInt(browser_version) > 8.0)){
        pdf_container =  "<iframe src='" + pdfjs_path + "' width='100%' height='750px'></iframe>";
        $('.highwire-pdf').html(pdf_container);
        $('.higwire-pdf-download-link').remove();
      }
      else if ((browser_name == "ie") (parseInt(browser_version) <=8.0)){
        show_pdf = true;
        if (acrobat) {
          $('.higwire-pdf-download-link').remove();
          var pdf_container = '<object width="100%"><embed src="' + pdf_path + '#view=fitbh&toolbars=1&pagemode=none&navpanes=0" width="100%" height="785px"></object>';
          $('.highwire-pdf').html(pdf_container);

          $(".ui-dialog").once("dialogopen", function(event, ui) {
            $(this).append('<iframe class="highwire-dialog-ie-wrap" src="about:blank"></iframe>');
            $('.highwire-dialog-ie-wrap').css('z-index','-1');
            $('.highwire-dialog-ie-wrap').css('position','absolute');
            $(this).css('z-index','100');
          });
        }
      }
      else {
        $('.highwire-pdf').remove();
      }
    }
  }
})(jQuery);

