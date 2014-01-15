
    var initialZoom = 0;
    var map;
    var centreLat = 0;
    var centreLon = -180;
	var markers = [];
	var mgr;
    var png_sizes = [];

    function pdfCoordToLatLng(x, y) {
        lng = (x / pdf_width - 0.5) * 360;
        lat = -((pdf_height - y) / pdf_height - 0.5) * 180;
        return new GLatLng(lat, lng);
    }
    function CustomProjection(zoom_level) {
      this.imageDimension = 65536;
      this.pixelsPerLonDegree = [];
      this.pixelsPerLatDegree = [];
      this.pixelOrigin = [];
      this.tileBounds = [];
      this.boundPixelSize = [];
      this.tileSize = 256;
      var cur_height = tile_height;
      var cur_width = tile_width;
      var tile_cntr = Math.ceil(Math.max(tile_height, tile_width) / 256);
      for(var d = 0; d < zoom_level; d++) {
        var mid_x = cur_width / 2;
        var mid_y = cur_height / 2;
        this.pixelsPerLonDegree.push(cur_width/360);
        this.pixelsPerLatDegree.push(cur_height/360);
        this.pixelOrigin.push(new GPoint(mid_x, mid_y));
        this.tileBounds.push(tile_cntr);
        png_sizes.push([cur_width, cur_height]);
        cur_width *= 2;
        cur_height *= 2;
        tile_cntr *= 2;
      }
    }
    CustomProjection.prototype = new GProjection();
    CustomProjection.prototype.fromLatLngToPixel = function(latlng, zoom) {
      var x = Math.round(
        this.pixelOrigin[zoom].x + 
          latlng.lng() * this.pixelsPerLonDegree[zoom]
      );
      var y = Math.round(
        this.pixelOrigin[zoom].y + 
          (-2 * latlng.lat()) * this.pixelsPerLatDegree[zoom]
      );
      return new GPoint(x, y);
    };
    CustomProjection.prototype.fromPixelToLatLng = function(pixel, zoom, unbounded) {
      var lng = (pixel.x - this.pixelOrigin[zoom].x) / this.pixelsPerLonDegree[zoom];
      var lat = -0.5*(pixel.y - this.pixelOrigin[zoom].y) / this.pixelsPerLatDegree[zoom];
      return new GLatLng(lat, lng, unbounded)
    };
    CustomProjection.prototype.tileCheckRange = function(tile, zoom, tilesize) {
      var tileBounds = this.tileBounds[zoom];
      if (tile.y < 0 || tile.y >= tileBounds)
        return false;
      if (tile.x < 0 || tile.x >= tileBounds)
        return false;
      return true;
    }
    CustomProjection.prototype.getWrapWidth = function(zoom) {
      return this.tileBounds[zoom] * this.tileSize;
    }

    function customGetTileURL(point, zoom_level) {
      var mid = Math.pow(2, zoom_level);
      var x = point.x;
      var y = point.y;
      return "/tree/" + pdffile + "/tile/?zoom="+zoom_level+"&tx=" + x + "&ty=" + y;
    }

    function load() {
      if (GBrowserIsCompatible()) {
        var copyright = new GCopyright(1,
                              new GLatLngBounds(new GLatLng(-90, -180),
                                                new GLatLng(90, 180)),
                              0,
                              "<a href=\"http://chemmine.ucr.edu\">ChemMine</a>");
        var copyrightCollection = new GCopyrightCollection("UCR IIGB");
        copyrightCollection.addCopyright(copyright);

        //create a custom picture layer
        var pic_tileLayers = [ new GTileLayer(copyrightCollection , 0, maxr)];
        pic_tileLayers[0].getTileUrl = customGetTileURL;
        pic_tileLayers[0].isPng = function() { return false; };
        pic_tileLayers[0].getOpacity = function() { return 1.0; };
        var proj = new CustomProjection(maxr+1);
        var pic_customMap = new GMapType(pic_tileLayers, proj, "Pic",
            {maxResolution:maxr, minResolution:0, errorMessage:"Data not available"});


        //Now create the custom map. Would normally be G_NORMAL_MAP,G_SATELLITE_MAP,G_HYBRID_MAP
        map = new GMap2(document.getElementById("map"),{mapTypes:[pic_customMap]});
        map.addControl(new GLargeMapControl());
        //map.addControl(new GMapTypeControl());
        map.addControl(new GOverviewMapControl(new GSize(300,300)));
        map.enableDoubleClickZoom();
        map.enableContinuousZoom();
        //map.enableScrollWheelZoom();
        map.setCenter(new GLatLng(centreLat, centreLon, true), initialZoom, pic_customMap);

		var _icon = new GIcon(G_DEFAULT_ICON);
		_icon.image = app_path + '/static/marker.png';
		_icon.shadow = app_path + '/static/shadow.png';
		_icon.shadowSize = new GSize(18, 18);
		_icon.iconSize = new GSize(18, 18);
		_icon.iconAnchor = new GPoint(18, 18);
		for (i in compounds) {
			v = compounds[i];
			m1 = new GMarker(pdfCoordToLatLng(v[0], v[1]), {icon:_icon});
			var compound = v[2]['label'];
			if (v[2]['url']) compound = "<a target=\"_blank\" href=\"" + v[2]['url'] + "\">" + v[2]['label'] + "</a>";
			if (typeof v[2]['starred'] == 'boolean') {
				if (v[2]['starred']) {
					compound = '<span class="star_ctrl"><a onclick="return star(this);" class="unstarbtn" href="' + app_path + '/work/compound/unstar/' + v[2]['label'] + '"><img src="' + app_path + '/static/star_on.png"/></a> <a onclick="return star(this)" class="starbtn" style="display:none" href="' + app_path + '/work/compound/star/' + v[2]['label'] + '"><img src="${app_path}/static/star_off.png"/></a></span>' + compound;
				} else {
					compound = '<span class="star_ctrl"><a onclick="return star(this)" class="unstarbtn" style="display:none" href="' + app_path + '/work/compound/unstar/' + v[2]['label'] + '"><img src="' + app_path + '/static/star_on.png"/></a> <a onclick="return star(this)" class="starbtn" href="' + app_path + '/work/compound/star/' + v[2]['label'] + '"><img src="' + app_path + '/static/star_off.png"/></a></span>' + compound;
				}
			}
			m1.bindInfoWindowHtml("<div style=\"min-width:250px; min-height:250px\">" + compound + "<br/><img src=\"" + app_path + "/work/compound/" + v[2]['label'] + "/png\"/></div>");
			markers.push(m1);
		}

		for (i in tree) {
			node = tree[i];
			if (!node[2]) continue;
			m1 = new GMarker(pdfCoordToLatLng(node[0], node[1]), {icon:_icon});
			m1.bindInfoWindowHtml("<div class=\"infowin\" style=\"min-width:250px; min-height:250px\">A branch of " + node[4]['n_leaves'] + " compounds.<br/><ul>" + 
			"<li><button onclick=\"branch(" + i + ")\"/>Show Compounds</button></li>" + 
			"<li><button onclick=\"branch(" + i + ", true)\"/>Show Compounds with structure</button></li>" + 
			"<li><button onclick=\"subtree('" + node[4]['repr_leaf'] + "'," + node[4]['height']+")\"/>Show only this branch</button></li>" + 
			"</ul>");
			markers.push(m1);
		}
		mgr = new MarkerManager(map);
		mgr.addMarkers(markers, marker_zoom, maxr);
		mgr.refresh();
		}

        GEvent.addListener(map, "zoomend", function(oz, nz) {
            jQuery("#downloadpng").text("Download " + png_sizes[nz][0] + "x" + png_sizes[nz][1] + " Image");
            jQuery("#downloadpng").siblings("label").text("You are at zoom level " + nz);
        });
    }

function branch(i, show_structure) {
	data = Array();
	function process_node(x) {
		x = tree[x];
		if (x[2] == null) {
			data.push(x[4]);
		} else {
			process_node(x[2]);
			process_node(x[3]);
		}
	}
	process_node(i);
	html = "<ul>";
	for (i in data) {
		html += "<li>"
		html += "<div>" + data[i]['label'] + "</div>";
		if (show_structure)
			html += "<div><img src=\"" + app_path + "/work/compound/" + data[i]['label'] + "/png\"/></div>";
		html += "</li>";
	}
	html += "</ul>";
	jQuery("#info").html(html).hide();
	jQuery("#info").css('max-height', window.innerHeight-230);
	jQuery("#infocontainer").show();
	jQuery("#info").slideDown();
	return false;
}

function show_wait()
{
    jQuery("#ctrl, #uploadform").hide();
	jQuery("#map").html("Rendering image. This can take up to a minute...");
}
function subtree(leaf, height)
{
    show_wait();
	document.location = document.location + '&subl=' + leaf + '&subh=' + height;
}

function star(me) {
	var me_ = jQuery(me);
	jQuery.getJSON(me_.attr('href'), function(data) {
		if (data['status'] == 'OK') {
			me_.parent().children("a").toggle();
		} else {
			alert('starring/unstarring failed!');
		}
	});
	return false;
}


var markeron = '';
jQuery(document).ready(function() {
	jQuery("#closeinfo").click(function() {
		jQuery("#info").slideUp();
		jQuery("#info").queue(function(){
			jQuery("#infocontainer").hide();
			jQuery(this).dequeue();
		});
		return false;
	});
	jQuery("#mlevel").text(marker_zoom);
	jQuery("#markerctrl button").click(function() {
		if (jQuery("#markerctrl").hasClass("on")) {
			mgr.clearMarkers();
			markeron = jQuery("#markerstatus").html();
			jQuery("#markerstatus").html("off");
			jQuery("#markerctrl button").text("enable markers");
			jQuery("#markerctrl").removeClass("on");
			jQuery("#markerctrl").addClass("off");
		} else {
			mgr.addMarkers(markers, marker_zoom, maxr);
			mgr.refresh();
			jQuery("#markerstatus").html(markeron);
			jQuery("#markerctrl button").text("disable markers");
			jQuery("#markerctrl").removeClass("off");
			jQuery("#markerctrl").addClass("on");

		}
	});
	jQuery(".showform").click(function() {
		jQuery("#uploadform").show("slow");
		return false;
	});
	jQuery("#closeform").click(function() {
		jQuery("#uploadform").hide("slow");
		return false;
	});
	jQuery("#treeonly").click(function() {
		loc = request_path + '?ref=' + ref;
        if (mode) loc += '&mode=' + mode;
        show_wait();
        document.location = loc;
	});
    jQuery("input[type=submit]").click(function() {
        show_wait();
        return true;
    });
    jQuery("#listcompounds").data('activated', false).click(function() {
        if (jQuery(this).data('activated')) {
            jQuery(this).data('activated', false);
            jQuery(this).text(jQuery(this).data('defaultlabel'));
            jQuery(this).siblings("textarea").remove();
        } else {
            var c = jQuery("<textarea readonly=\"readonly\" style=\"display:block; width:100%; height:150px; overflow:auto\"></textarea>").insertAfter(jQuery(this));
            var content = '';
            for (i in compounds) {
                content += compounds[i][2]['label'] + "\n";
            }
            c.text(content).focus().select();
            c.focus(function() {jQuery(this).select();});
            jQuery(this).data('defaultlabel', jQuery(this).text());
            jQuery(this).text('Hide List');
            jQuery(this).data('activated', true);
        }
    }).before("<label>" + compounds.length + " compounds appear in this " + (mode=='sp'?"plot":"tree") + "</label>");
    jQuery("#ctrlctrl").data('minimized', false);
    jQuery("#ctrlctrl").data('boxwidth', jQuery("#ctrl").get(0).offsetWidth);
    jQuery("#ctrlctrl").data('boxheight', jQuery("#ctrl").get(0).offsetHeight);
    jQuery("#ctrlctrl").attr('title', 'minimize control box');
    jQuery("#ctrlctrl").click(function() {
        if (! jQuery(this).data('minimized')) {
            jQuery(this).children("span:eq(0)").hide();
            jQuery(this).children("span:eq(1)").show();
            jQuery(this).parent().queue(function() {
                jQuery(this).children("div, hr").hide(100);
                jQuery(this).dequeue();
            });
            jQuery(this).parent().animate({width:'14px', height:'14px'}, 100);
            jQuery(this).data('minimized', true);
            jQuery(this).attr('title', 'show control box');
        }
        else {
            jQuery(this).children("span:eq(1)").hide();
            jQuery(this).children("span:eq(0)").show();
            jQuery(this).parent().animate({width:jQuery(this).data('boxwidth'), height:jQuery(this).data('boxheight')}, 100, null, function() {
                jQuery(this).children("div, hr").show();
            });
            jQuery(this).data('minimized', false);
            jQuery(this).attr('title', 'minimize control box');
        }
        return false;
    });
    jQuery("#downloadpng").click(function() {
        if (confirm("This will open an image of size " + png_sizes[map.getZoom()][0] + "x" + png_sizes[map.getZoom()][1] + " in your browser. Continue?")) {
            document.location = app_path + '/tree/'+ pdffile + '/full?zoom=' + map.getZoom();
        }
    });
    jQuery("#downloadpng").text("Download " + png_sizes[0][0] + "x" + png_sizes[0][1] + " Image").before("<label>You are at zoom level " + map.getZoom() + "</label>");
});

