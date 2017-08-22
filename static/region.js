var genome_coords_margin = {left: 80, right: 30};
var coverage_plot_height = 200;
var coverage_plot_margin = {top: 10, bottom: 30};
var variant_plot_height = 30;
var variant_plot_margin = {top: 5, bottom: 5};

function bootstrap_plot() {
    window.model = window.model || {};
    window.model.plot = window.model.plot || {};
    if (!window.model.plot.svg_width) {
        window.model.plot.svg_width = $('#coverage_plot_container').width();
        window.model.plot.genome_coords_width = window.model.plot.svg_width - genome_coords_margin.left - genome_coords_margin.right;
    }
    if (!window.model.plot.x) {
        window.model.plot.x = d3.scale.linear()
            .domain([window.model.start, window.model.stop])
            .range([0, window.model.plot.genome_coords_width]);
    }
}

function create_coverage_chart(cov_data) {
    bootstrap_plot();

    var metric = 'mean';

    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(cov_data, function(d) { return d[metric]; });
    }

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([coverage_plot_height, 0]);

    var xAxis = d3.svg.axis()
        .scale(window.model.plot.x)
        .orient("bottom");

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select('#coverage_plot_container').append("svg")
        .attr("width", window.model.plot.svg_width)
        .attr("height", coverage_plot_height + coverage_plot_margin.top + coverage_plot_margin.bottom)
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + genome_coords_margin.left + "," + coverage_plot_margin.top + ")");

    svg.selectAll("bar")
        .data(cov_data)
        .enter()
        .append("rect")
        .attr('class', 'main_plot_bars')
        .style("fill", "steelblue")
        .attr("x", function(d) {
            return window.model.plot.x(d.start);
        })
        .attr("width", function(d) {
            var length_in_bases = d.end - d.start + 1;
            var width_of_base = window.model.plot.genome_coords_width/cov_data.length;
            return length_in_bases * width_of_base;
        })
        .attr("y", function(d) {
            return y(d[metric]) || 0;
        })
        .attr("height", function(d) {
            return (coverage_plot_height - y(d[metric])) || 0;
        });

    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + coverage_plot_height + ")")
        .call(xAxis);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);
}

function change_coverage_chart_metric(cov_data, metric) {

    console.log(['cov_data', cov_data]);

    var max_cov = 1;
    if (metric === 'mean' || metric === 'median') {
        max_cov = d3.max(cov_data, function(d) { return d[metric]; });
    }
    console.log(max_cov);

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([coverage_plot_height, 0]);

    var svg = d3.select('#coverage_plot_container').select('#inner_graph');

    svg.selectAll("rect")
        .data(cov_data)
        .transition()
        .duration(500)
        .attr("y", function(d) { return y(d[metric]); })
        .attr("height", function(d) { return coverage_plot_height - y(d[metric]); });

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);
}

function create_variant_plot() {
    bootstrap_plot();

    var svg = d3.select('#variant_plot_container').append("svg")
        .attr("width", window.model.plot.svg_width)
        .attr("height", variant_plot_height)
        .append("g")
        .attr('id', 'variant_track')
        .attr("transform", "translate(" + genome_coords_margin.left + "," + 0 + ")");

    svg.append("line")
        .attr("y1", variant_plot_height/2)
        .attr("y2", variant_plot_height/2)
        .attr("x1", 0)
        .attr("x2", window.model.plot.genome_coords_width)
        .attr("stroke-width", 1)
        .attr("stroke", "lightsteelblue");

    window.model.plot.oval_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        var csq = d.major_consequence.replace(/_/g, ' ');
        if (csq.length > 15) { csq = csq.substr(0, 15) + '...'; } // because d3-tip tooltips fall off the page
        return group_thousands(d.pos) + '<br>' +
            csq + '<br>' +
            (d.filter === 'PASS' ? '' : d.filter + '<br>') +
            'MAF: ' + perc_sigfigs(d.allele_freq, 2);
        //return JSON.stringify(d);
    });
    svg.call(window.model.plot.oval_tip);
}

function change_variant_plot(variants) {
    //console.log('plotting', variants);
    var selection = d3.select('#variant_track').selectAll('.variant-circle').data(variants);
    selection.enter() //ENTER
        .append('ellipse')
        //.attr('foo', function(d) { console.log('enter-each', d)})
        .attr('class', 'variant-circle')
        .style('opacity', 0.3)
        .style('fill', 'blue')
        .attr('ry', 4)
        .attr('rx', 4)
        .attr('cy', variant_plot_height/2)
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
        .on('mouseover', window.model.plot.oval_tip.show)
        .on('mouseout', window.model.plot.oval_tip.hide)
    selection // UPDATE // todo: learn how to get union of ENTER+UPDATE
        //.attr('foo', function(d) { console.log('update-each', d)})
        .attr('cx', function(d) { return window.model.plot.x(d.pos); })
    selection.exit() //EXIT
        //.attr('foo', function(d) { console.log('exit-each', d)})
        .remove()
}


$(document).ready(function() {
    if (window.model.coverage_stats != null) {
        create_coverage_chart(window.model.coverage_stats);
        // Change coverage plot
        $('.coverage_metric_buttons').change(function () {
            var v = $(this).attr('id').replace('_covmet_button', '');
            $('.coverage_subcat_selectors').hide();
            if (v == 'covered') {
                $('#over_x_select_container').show();
                v = $('#over_x_select').val().replace('X', '');
            } else {
                $('#average_select_container').show();
                v = $("#average_select").val();
            }
            change_coverage_chart_metric(window.model.coverage_stats, v);
        });
        $('#over_x_select').change(function () {
            change_coverage_chart_metric(window.model.coverage_stats, $(this).val().replace('X', ''));
        });
        $('#average_select').change(function () {
            change_coverage_chart_metric(window.model.coverage_stats, $(this).val());
        });
    }

    create_variant_plot();
});


function create_variant_table() {
    var columns = [
        {
            title: 'Alleles (rsID)', name: 'allele',
            searchable: false, orderable: false,
            render: (function() {
                var template = _.template(
                    '<a href="/variant/<%= window.model.chrom %>-<%= variant.pos %>-<%= variant.ref %>-<%= variant.alt %>" target="_blank">'+
                        '<% if (variant.ref.length > 20) { %><%= variant.ref.slice(0,20) %>...<% } else { %><%= variant.ref %><% } %> / '+
                        '<% if (variant.alt.length > 20) { %><%= variant.alt.slice(0,20) %>...<% } else { %><%= variant.alt %><% } %>'+
                        '<% if (variant.rsids.length) { %>(<%= variant.rsids.join(", ") %>)<% } %>'+
                        '</a>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Position on chr' + window.model.chrom, name: 'pos',
            data: 'pos', searchable: true,  orderable: true, className: 'dt-right',
            render: function(cell_data, type, row) { return group_thousands(cell_data); },

        },{
            title: 'HGVS', name: 'hgvs',
            searchable: false, orderable: false, className: 'dt-head-center',
            render: (function() {
                var template = _.template(
                    '<% if (variant.HGVS != "") { %>' +
                    '<% if (variant.HGVS.length > 20) { %><%= variant.HGVS.slice(0,20) %>...<% } else { %><%= variant.HGVS %><% } %>' +
                        '<% if (variant.CANONICAL !== "YES") { %><span class="tooltip-table-header" data-tooltip="Annotation is for non-canonical transcript">&dagger;</span><% } %>'+
                        '<% } %>',
                    {variable:'variant'});
                return function(cell_data, type, row) { return template(row); };
            })(),

        },{
            title: 'Annotation', name: 'csq',
            data: 'major_consequence', searchable:true, orderable:false, className: 'dt-center',
            render: function(cell_data, type, row) {
                return fmt_annotation(cell_data);
            },

        },{
            title: 'CADD', name:'cadd',
            orderable:false, className: 'dt-right',
            render: function() { return '0'; },

        },{
            title: 'QC', name: 'filter',
            data: 'filter', searchable:true, orderable:false, className: 'dt-center',

        },{
            title: 'HomRef', name: 'n_homref',
            searchable:true, orderable:false, className: 'dt-right',
            render: function(cell_data, type, row) {
                var num_het_samples = row.allele_count - 2*row.hom_count;
                var num_homref_samples = ((row.allele_num - row.allele_count) - num_het_samples) / 2;
                return group_thousands(num_homref_samples);
            },

        },{
            title: 'Het', name: 'n_het',
            searchable:true, orderable:false, className: 'dt-right',
            render: function(cell_data, type, row) {
                var num_het_samples = row.allele_count - 2*row.hom_count;
                return group_thousands(num_het_samples);
            },

        },{
            title: 'HomAlt', name: 'hom_count',
            searchable:true, orderable:true, className: 'dt-right',
            render: function(cell_data, type, row) {
                return group_thousands(row.hom_count);
            },

        },{
            title: 'Frequency', name: 'allele_freq',
            data: 'allele_freq', searchable:true, orderable:true, className: 'dt-pad-left',
            render: function(cell_data, type, row) { return perc_sigfigs(cell_data, 2); },

        },
    ];

    window.model = window.model || {};
    window.model.filter_info = window.model.filter_info || {};
    window.model.filter_info.start = window.model.start;
    window.model.filter_info.stop = window.model.stop;
    window.model.filter_info.chrom = window.model.chrom;

    $('#variant_table').DataTable({
        serverSide: true, /* API does all the real work */

        processing: true, /* show "processing" over table while waiting for API */
        deferRender: true, /* only render rows when they're being displayed */

        paging: true,
        pagingType: 'full', /* [first, prev, next, last] */
        pageLength: 10,

        searching: false,

        ajax: {
            url: '/api/variants_for_table',
            type: 'POST',
            data: function(args) { /* modify request form parameters */
                return {
                    args: JSON.stringify(args), // jsonify all params rather than using `columns[0][search][value]` php form syntax
                    filter_info: JSON.stringify(window.model.filter_info),
                };
            },
            dataSrc: function(resp) { /* modify API's response (for debugging) */
                window._debug = window._debug || {}; window._debug.resp = resp;
                change_variant_plot(resp.data);
                return resp.data;
            },
        },

        order: [[columns.map(function(d){return d.data==='pos'}).indexOf(true), 'asc']],
        columns: columns,

        dom: '<ip>ftr'  //'ipftr',  // default is 'lfrtip'.  l=length f=filtering t=table i=info p=paging, r=processing

    });

    $('.variant_table_filter').on('change', function() {
        window.model.filter_info.pos_ge = parseInt($('input#pos_ge').val());
        window.model.filter_info.pos_le = parseInt($('input#pos_le').val());
        window.model.filter_info.filter_value = $('select#filter_value').val();
        console.log(window.model);
        $('#variant_table').DataTable().draw();
    });
}
$(create_variant_table);