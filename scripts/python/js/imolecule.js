/*global THREE, $, jQuery, console, window, requestAnimationFrame, document, io, alert, Blob, saveAs*/
"use strict";

var imolecule = {

    // Creates a new instance of imolecule
    create: function (selector, options) {
        var $s = $(selector), self = this;
        options = options || {};

        this.shader = options.hasOwnProperty("shader") ? options.shader : THREE.ShaderToon.toon2;
        this.drawingType = options.hasOwnProperty("drawingType") ? options.drawingType : "ball and stick";
        this.boundaryType = options.hasOwnProperty("boundaryType") ? options.boundaryType : "unit cell";
        this.renderer = new THREE.WebGLRenderer({antialias: true});
        this.renderer.setSize($s.width(), $s.height());
        $s.append(this.renderer.domElement);

        this.perspective = new THREE.PerspectiveCamera(50, $s.width() / $s.height());
        this.orthographic = new THREE.OrthographicCamera(-$s.width() / 32,
                $s.width() / 32, $s.height() / 32, -$s.height() / 32, -100, 1000);
        this.perspective.position.z = options.hasOwnProperty("z") ? options.z : 15;
        this.orthographic.position = this.perspective.position;
        this.orthographic.rotation = this.perspective.rotation;
        this.setCameraType(options.hasOwnProperty("cameraType") ? options.cameraType : "perspective");

        this.sphereGeometry = new THREE.SphereGeometry(1, 16, 12);
        this.cylinderGeometry = new THREE.CylinderGeometry(1, 1, 1, 6, 3, false);

        // This orients the cylinder primitive so THREE.lookAt() works properly
        this.cylinderGeometry.applyMatrix(new THREE.Matrix4()
            .makeRotationFromEuler(new THREE.Vector3(Math.PI / 2, Math.PI, 0)));

        this.light = new THREE.HemisphereLight(0xffffff, 1.0);
        this.light.position = this.camera.position;
        this.light.rotation = this.camera.rotation;

        this.atoms = [];
        this.bonds = [];

        // Makes toon-shaded materials from scratch
        $.each(this.data, function (key, value) {
            var material = new THREE.ShaderMaterial({
                uniforms: THREE.UniformsUtils.clone(self.shader.uniforms),
                vertexShader: self.shader.vertexShader,
                fragmentShader: self.shader.fragmentShader
            });
            material.uniforms.uDirLightPos.value.set(1, 1, 1)
                             .multiplyScalar(self.camera.position.z);
            value.color = new THREE.Color(value.color);
            material.uniforms.uDirLightColor.value = value.color;
            material.uniforms.uBaseColor.value = value.color;
            value.material = material;
        });

        // Initializes a scene and appends objects to be drawn
        this.scene = new THREE.Scene();
        this.scene.add(this.perspective);
        this.scene.add(this.orthographic);
        this.scene.add(this.light);

        $(window).resize(function () {
            self.renderer.setSize($s.width(), $s.height());
            self.perspective.aspect = $s.width() / $s.height();
            self.perspective.updateProjectionMatrix();
            self.orthographic.left = -$s.width() / 32.0;
            self.orthographic.right = $s.width() / 32.0;
            self.orthographic.top = $s.height() / 32.0;
            self.orthographic.bottom = -$s.height() / 32.0;
            self.orthographic.updateProjectionMatrix();
        });

        this.animate();
    },

    // Draws a molecule. Duh.
    draw: function (molecule) {
        var mesh, self, a, scale, j, k, dy, cent, data, v, vectors, points,
            trans, geometry, material;
        self = this;
        cent = new THREE.Vector3();
        this.current = molecule;

        scale = this.drawingType === "space filling" ? 1.0 : 0.3;

        // Don't hate on formats without bond information
        if (!molecule.hasOwnProperty("bonds")) { molecule.bonds = []; }

        // Draws atoms and saves references
        $.each(molecule.atoms, function (i, atom) {
            data = self.data[atom.element] || self.data.unknown;
            mesh = new THREE.Mesh(self.sphereGeometry, data.material);
            mesh.position.fromArray(atom.location);
            mesh.scale.set(1, 1, 1).multiplyScalar(scale * data.radius * 2);
            if (self.drawingType !== "wireframe") {
                self.scene.add(mesh);
            }
            mesh.element = atom.element;
            self.atoms.push(mesh);
        });

        // Bonds require some basic vector math
        $.each(molecule.bonds, function (i, bond) {
            a = [self.atoms[bond.atoms[0]], self.atoms[bond.atoms[1]]];
            for (j = 0; j < bond.order; j += 1) {
                if (bond.order === 2) {
                    dy = 0.5 * ((j === 1) ? 1 : -1);
                } else if (bond.order === 3 && j !== 0) {
                    dy = ((j === 1) ? 1 : -1);
                } else {
                    dy = 0;
                }

                for (k = 0; k < 2; k += 1) {
                    mesh = new THREE.Mesh(self.cylinderGeometry, self.data.bond.material);
                    cent.addVectors(a[0].position, a[1].position).divideScalar(2);
                    mesh.atomMaterial = self.data[a[k].element].material;
                    mesh.position.addVectors(cent, a[k].position).divideScalar(2);
                    mesh.lookAt(a[1].position);
                    mesh.scale.x = mesh.scale.y = 0.3 * self.data.bond.radius * 2;
                    mesh.scale.z = a[1].position.distanceTo(a[0].position) / 2.0;
                    mesh.translateY(0.3 * dy);

                    if (self.drawingType === "wireframe") {
                        mesh.material = mesh.atomMaterial;
                    }
                    if (self.drawingType !== "space filling") {
                        self.scene.add(mesh);
                    }
                    self.bonds.push(mesh);
                }
            }
        });

        // If we're dealing with a crystal structure, draw the unit cell
        if (molecule.hasOwnProperty("periodic_connections")) {
            // Some basic conversions to handle math via THREE.Vector3
            v = new THREE.Vector3(0, 0, 0);
            vectors = [
                v.clone().fromArray(molecule.periodic_connections[0]),
                v.clone().fromArray(molecule.periodic_connections[1]),
                v.clone().fromArray(molecule.periodic_connections[2])
            ];
            // The eight corners of the unit cell are linear combinations of above
            points = [
                v.clone(), vectors[0], vectors[1], vectors[2],
                v.clone().add(vectors[0]).add(vectors[1]).add(vectors[2]),
                v.clone().add(vectors[1]).add(vectors[2]),
                v.clone().add(vectors[0]).add(vectors[2]),
                v.clone().add(vectors[0]).add(vectors[1])
            ];
            // Translate unit cell to center around mof + origin
            trans = points[4].clone().multiplyScalar(0.5);
            for (j = 0; j < points.length; j += 1) {
                points[j].sub(trans);
            }
            // Draw the box line-by-line
            geometry = new THREE.Geometry();
            $.each([0, 1, 0, 2, 0, 3, 6, 1, 7, 2, 5, 3, 5, 4, 6, 4, 7], function (index, value) {
                geometry.vertices.push(points[value]);
            });
            material = new THREE.LineBasicMaterial({color: 0x000000, linewidth: 3});
            this.corners = new THREE.Line(geometry, material);
            this.scene.add(this.corners);
        }
    },

    // Deletes any existing molecules.
    clear: function () {
        var self = this;
        $.each(this.atoms.concat(this.bonds), function (i, value) {
            self.scene.remove(value);
        });
        this.atoms = [];
        this.bonds = [];
        this.scene.remove(this.corners);
    },

    // Sets molecule drawing types ( ball and stick, space filling, wireframe )
    setDrawingType: function (type) {
        // Some case-by-case logic to avoid clearing and redrawing the canvas
        var i;
        if (this.drawingType === "ball and stick") {
            if (type === "wireframe") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.scene.remove(this.atoms[i]);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.bonds[i].material = this.bonds[i].atomMaterial;
                }
            } else if (type === "space filling") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.atoms[i].scale.divideScalar(0.3);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.scene.remove(this.bonds[i]);
                }
            }
        } else if (this.drawingType === "wireframe") {
            if (type === "ball and stick") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.scene.add(this.atoms[i]);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.bonds[i].material = this.data.bond.material;
                }
            } else if (type === "space filling") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.atoms[i].scale.divideScalar(0.3);
                    this.scene.add(this.atoms[i]);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.scene.remove(this.bonds[i]);
                }
            }
        } else if (this.drawingType === "space filling") {
            if (type === "ball and stick") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.atoms[i].scale.multiplyScalar(0.3);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.bonds[i].material = this.data.bond.material;
                    this.scene.add(this.bonds[i]);
                }
            } else if (type === "wireframe") {
                for (i = 0; i < this.atoms.length; i += 1) {
                    this.atoms[i].scale.multiplyScalar(0.3);
                    this.scene.remove(this.atoms[i]);
                }
                for (i = 0; i < this.bonds.length; i += 1) {
                    this.bonds[i].material = this.bonds[i].atomMaterial;
                    this.scene.add(this.bonds[i]);
                }
            }
        }
        this.drawingType = type;
    },

    // Sets camera type (orthogonal, perspective)
    setCameraType: function (type) {
        if (type === "orthographic") {
            this.camera = this.orthographic;
            this.camera.position.copy(this.perspective.position);
        } else if (type === "perspective") {
            this.camera = this.perspective;
            this.camera.position.copy(this.orthographic.position);
        }
        this.controls = new THREE.TrackballControls(this.camera, this.renderer.domElement);
    },

    // Runs the main window animation in an infinite loop
    animate: function () {
        var self = this;
        window.requestAnimationFrame(function () {
            return self.animate();
        });
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    },

    // Either shows or hides the unit cell
    showUnitCell: function (toggle) {
        this.scene[toggle ? "add" : "remove"](this.corners);
    },

    // Connects to Python via a socketio-zeromq bridge. Ignore everything below
    // if you're not using the client-server functionality
    connect: function (http_port) {
        this.socket = new io.Socket(window.location.hostname,
                                    {port: http_port, rememberTransport: false});
        this.socket.connect();
        this.queue = {};

        this.socket.on("connect", function () {
            console.log("Connected!");
        });

        var self = this;
        this.socket.on("message", function (data) {
            var router, name;
            router = self.queue[data.id];
            delete self.queue[data.id];
            self.result = data.result;
            try { data.result = $.parseJSON(data.result); } catch (err) {}

            if (data.error) {
                alert(data.result);
                return;
            }

            if (router === "draw") {
                self.clear();
                self.draw(data.result);
            } else if (router === "save") {
                name = (data.result.hasOwnProperty("name") && data.result.name !== "") ?
                        data.result.name : self.filename;
                saveAs(new Blob([self.result], {type: "text/plain"}),
                        name + "." + self.outFormat);
            } else {
                alert("Unsupported function: " + router);
            }
        });
    },

    // Standardizes chemical drawing formats and draws
    convertAndDraw: function (data, inFormat) {
        // Skips a socket call if not needed
        if (inFormat === "json") {
            this.clear();
            this.draw($.parseJSON(data));
            return;
        }

        var uuid = this.uuid();
        this.socket.send({method: "convert", id: uuid,
                          params: {data: data, in_format: inFormat,
                                   out_format: "json", pretty: false}
                         });
        this.queue[uuid] = "draw";
    },

    // Converts and saves output
    convertAndSave: function (data, outFormat) {
        var uuid = this.uuid();
        this.socket.send({method: "convert", id: uuid,
                          params: {data: data, in_format: "json",
                                   out_format: outFormat, pretty: true}
                         });
        this.outFormat = outFormat;
        this.queue[uuid] = "save";
    },

    // Generates a unique identifier for request ids
    // Code from http://stackoverflow.com/questions/105034/
    // how-to-create-a-guid-uuid-in-javascript/2117523#2117523
    uuid: function () {
        return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, function (c) {
            var r = Math.random()*16|0, v = c == 'x' ? r : (r&0x3|0x8);
            return v.toString(16);
        });
    },

    data: { Ac: { color: 0x70aaf9, radius: 1.95 },
            Ag: { color: 0xbfbfbf, radius: 1.6 },
            Al: { color: 0xbfa5a5, radius: 1.25 },
            Am: { color: 0x545bf2, radius: 1.75 },
            Ar: { color: 0x80d1e2, radius: 0.71 },
            As: { color: 0xbc80e2, radius: 1.15 },
            Au: { color: 0xffd123, radius: 1.35 },
            B: { color: 0xffb5b5, radius: 0.85 },
            Ba: { color: 0x00c800, radius: 2.15 },
            Be: { color: 0xc1ff00, radius: 1.05 },
            Bi: { color: 0x9e4fb5, radius: 1.6 },
            Br: { color: 0xa52828, radius: 1.15 },
            C: { color: 0x909090, radius: 0.7 },
            Ca: { color: 0x3dff00, radius: 1.8 },
            Cd: { color: 0xffd88e, radius: 1.55 },
            Ce: { color: 0xffffc6, radius: 1.85 },
            Cl: { color: 0x1fef1f, radius: 1.0 },
            Co: { color: 0xef90a0, radius: 1.35 },
            Cr: { color: 0x8999c6, radius: 1.4 },
            Cs: { color: 0x56178e, radius: 2.6 },
            Cu: { color: 0xc88033, radius: 1.35 },
            Dy: { color: 0x1fffc6, radius: 1.75 },
            Er: { color: 0x00e675, radius: 1.75 },
            Eu: { color: 0x60ffc6, radius: 1.85 },
            F: { color: 0x90df4f, radius: 0.5 },
            Fe: { color: 0xdf6633, radius: 1.4 },
            Ga: { color: 0xc18e8e, radius: 1.3 },
            Gd: { color: 0x44ffc6, radius: 1.8 },
            Ge: { color: 0x668e8e, radius: 1.25 },
            H: { color: 0xffffff, radius: 0.25 },
            Hf: { color: 0x4dc1ff, radius: 1.55 },
            Hg: { color: 0xb8b8cf, radius: 1.5 },
            Ho: { color: 0x00ff9c, radius: 1.75 },
            I: { color: 0x930093, radius: 1.4 },
            In: { color: 0xa57572, radius: 1.55 },
            Ir: { color: 0x175487, radius: 1.35 },
            K: { color: 0x8e3fd4, radius: 2.2 },
            La: { color: 0x70d4ff, radius: 1.95 },
            Li: { color: 0xcc80ff, radius: 1.45 },
            Lu: { color: 0x00aa23, radius: 1.75 },
            Mg: { color: 0x89ff00, radius: 1.5 },
            Mn: { color: 0x9c79c6, radius: 1.4 },
            Mo: { color: 0x54b5b5, radius: 1.45 },
            N: { color: 0x2f4ff7, radius: 0.65 },
            Na: { color: 0xaa5bf2, radius: 1.8 },
            Nb: { color: 0x72c1c8, radius: 1.45 },
            Nd: { color: 0xc6ffc6, radius: 1.85 },
            Ni: { color: 0x4fcf4f, radius: 1.35 },
            Np: { color: 0x0080ff, radius: 1.75 },
            O: { color: 0xff0d0d, radius: 0.6 },
            Os: { color: 0x266695, radius: 1.3 },
            P: { color: 0xff8000, radius: 1.0 },
            Pa: { color: 0x00a1ff, radius: 1.8 },
            Pb: { color: 0x565960, radius: 1.8 },
            Pd: { color: 0x006985, radius: 1.4 },
            Pm: { color: 0xa3ffc6, radius: 1.85 },
            Po: { color: 0xaa5b00, radius: 1.9 },
            Pr: { color: 0xd8ffc6, radius: 1.85 },
            Pt: { color: 0xcfcfdf, radius: 1.35 },
            Pu: { color: 0x006bff, radius: 1.75 },
            Ra: { color: 0x007c00, radius: 2.15 },
            Rb: { color: 0x702daf, radius: 2.35 },
            Re: { color: 0x267caa, radius: 1.35 },
            Rh: { color: 0x0a7c8c, radius: 1.35 },
            Ru: { color: 0x238e8e, radius: 1.3 },
            S: { color: 0xffff2f, radius: 1.0 },
            Sb: { color: 0x9e62b5, radius: 1.45 },
            Sc: { color: 0xe6e6e6, radius: 1.6 },
            Se: { color: 0xffa100, radius: 1.15 },
            Si: { color: 0xefc8a0, radius: 1.1 },
            Sm: { color: 0x8effc6, radius: 1.85 },
            Sn: { color: 0x668080, radius: 1.45 },
            Sr: { color: 0x00ff00, radius: 2.0 },
            Ta: { color: 0x4da5ff, radius: 1.45 },
            Tb: { color: 0x2fffc6, radius: 1.75 },
            Tc: { color: 0x3b9e9e, radius: 1.35 },
            Te: { color: 0xd47900, radius: 1.4 },
            Th: { color: 0x00baff, radius: 1.8 },
            Ti: { color: 0xbfc1c6, radius: 1.4 },
            Tl: { color: 0xa5544d, radius: 1.9 },
            Tm: { color: 0x00d452, radius: 1.75 },
            U: { color: 0x008eff, radius: 1.75 },
            V: { color: 0xa5a5aa, radius: 1.35 },
            W: { color: 0x2193d6, radius: 1.35 },
            Y: { color: 0x93ffff, radius: 1.8 },
            Yb: { color: 0x00bf38, radius: 1.75 },
            Zn: { color: 0x7c80af, radius: 1.35 },
            Zr: { color: 0x93dfdf, radius: 1.55 },
            bond: { color: 0x0c0c0c, radius: 0.18 },
            unknown: { color: 0x000000, radius: 0.8 }
        }
};
