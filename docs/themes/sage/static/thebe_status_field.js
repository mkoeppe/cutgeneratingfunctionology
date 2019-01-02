// Functions for the thebe activate button and status field
//
// 

function thebe_place_activate_button(){
  $('.thebe_status_field')
    .html('<input type="button"\
                  onclick="thebe_activate_button_function()"\
                  value="Activate"\
                  title="ThebeLab (requires internet):\nClick to activate code cells in this page.\nYou can then edit and run them.\nComputations courtesy of mybinder.org."\
                   class="thebe-status-field"/>');
}

function thebe_remove_activate_button(){
  $('.thebe_status_field').empty();
}

function thebe_place_status_field(){
  $('.thebe_status_field')
    .html('<div class="thebe-status-field"\
                title="ThebeLab status.\nClick `run` to execute code cells.\nComputations courtesy of mybinder.org.">\
          </div>');
}
sagemath='<a href="http://sagemath.org">Sage</a>'
mybinder='<a href="http://mybinder.org">mybinder.org</a>'
about=' (<a href="http://sage-package.readthedocs.io/en/latest/sage_package/thebe.html">About</a>)'
server=''

messages = {
  'downloading': 'Downloading Thebe',
  'building': 'Building '+sagemath,
  'built': 'Built '+sagemath,
  'launching': 'Launching server',
  'server-ready': 'Launched server',
  'starting': 'Launching '+sagemath,
  'ready': 'Running  '+sagemath,
  'failed': 'Failed (network error?)'
}

function thebe_update_status_field(evt, data) {
  console.log("Thebe: status changed (" + data.status + "): " + data.message);
  $(".thebe-status-field")
    .attr("class", "thebe-status-field thebe-status-" + data.status)
    .html(messages[data.status]+server+about);
}

function thebe_bootstrap_local () {
  console.log("Thebe: using local server");
  thebelab.on("status", thebe_update_status_field);
  thebelab.bootstrap({
    binderOptions: {repo: null},
    kernelOptions: {
      name: "sagemath",
      serverSettings: {
        "baseUrl": window.location.origin,
      }
    }
  });
}

function thebe_bootstrap_binder () {
  console.log("Thebe: using remote server on binder");
  thebelab.on("status", thebe_update_status_field); // Duplicated with above; would be nicer as thebe option
  server = " on "+mybinder;
  thebelab.bootstrap();
}

// Try downloading thebe remotely; if successfull call next_operation
function thebe_download_remote(next_operation) {
  thebe_update_status_field({}, {status: 'downloading', message: ''})
  // Load the Thebe library
  $.getScript("https://unpkg.com/thebelab@^0.1.0")
    .done(function(script, textStatus) {
      next_operation()
    })
    .fail(function(jqxhr, settings, exception ) {
      $( "div.log" ).text( "Could not fetch ThebeLab library." );
      thebe_update_status_field({}, {status: 'failed', message: ''})
    });
}

// Try downloading thebe locally, or remotely if unavailable; if successfull call next_operation
function thebe_download_local(next_operation) {
  console.log("Thebe: trying to get thebe from the nbextensions");
  $.getScript("/nbextensions/thebelab.js")
    .done(function(script, textStatus) {
      next_operation()
    })
    .fail(function(jqxhr, settings, exception ) {
      thebe_download_remote(next_operation);
    });
}

// Activate button function hook
function thebe_activate_button_function(){
  thebe_remove_activate_button();
  thebe_place_status_field();
  // Checks whether served by a Jupyter server, and whether a sagemath kernel is available
  // Proceed accordingly by downloading thebe locally or remotely
  // and running a kernel locally or remotely
  if (window.location.protocol.startsWith('http')) {
    ajax(window.location.origin+'/api/kernelspecs', {
      dataType: 'json',
      success: function(json) {
        if (json['kernelspecs']['sagemath'])
          thebe_download_local(thebe_bootstrap_local)
      },
      error:    function() {
        thebe_download_local(thebe_bootstrap_binder);
      }
    });
  } else {
    thebe_download_remote(thebe_bootstrap_binder);
  }
}

/*****************************************************************************
  Jupyterlab utilities, taken from:
  https://github.com/jupyter/notebook/blob/master/notebook/static/base/js/utils.js

   Might not be needed once jupyterlab will be shiped with SageMath.
*****************************************************************************/

var _get_cookie = function (name) {
  // from tornado docs: http://www.tornadoweb.org/en/stable/guide/security.html
  var r = document.cookie.match("\\b" + name + "=([^;]*)\\b");
  return r ? r[1] : undefined;
}

var _add_auth_header = function (settings) {
  /**
   * Adds auth header to jquery ajax settings
   */
  settings = settings || {};
  if (!settings.headers) {
    settings.headers = {};
  }
  if (!settings.headers.Authorization) {
    var xsrf_token = _get_cookie('_xsrf');
    if (xsrf_token) {
      settings.headers['X-XSRFToken'] = xsrf_token;
    }
  }
  return settings;
};

var ajax = function (url, settings) {
  // like $.ajax, but ensure XSRF or Authorization header is set
  if (typeof url === "object") {
    // called with single argument: $.ajax({url: '...'})
    settings = url;
    url = settings.url;
    delete settings.url;
  }
  settings = _add_auth_header(settings);
  return $.ajax(url, settings);
};

