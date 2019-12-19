// DGP 2019 Project
// ShapeUp and Projective Dynamics
// Author: Christopher Brandt

// An API that provides simple methods to handle the interaction between
// the viewer and the ProjectiveDynamics simulation

#pragma once

#include "projdyn_types.h"
#include <thread>
#include "projdyn.h"
#include "projdyn_tetgen.h"
#include <nanogui/slider.h>
#include <nanogui/checkbox.h>
#include "viewer.h"
#include "projdyn_widgets.h"

#include <surface_mesh/types.h>
#include <surface_mesh/Surface_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

class ProjDynAPI {
public:
    // Some default values
    const bool UPDATE_NORMALS = true;
    const int  NUM_ITS_INITIAL = 10;
    const int  FPS = 60;
    const bool DYNAMIC_MODE = true;

    /** NEW values for diffusion : button toggling */
    bool diff_activated = false;
    enum TEMP_DIFFUSION_TYPE : int { UNILAPLACE = 2, COTANLAPLACE = 3};
    TEMP_DIFFUSION_TYPE diffusion_type = UNILAPLACE;
    float decay = 0.0f;
    int diffusion_iterations = 1;
    ///////////////////////////////////

    Eigen::VectorXd m_temperatures;

    bool init_temp = false;

    bool is_tetra =false;


    ProjDynAPI(Viewer* viewer) {
        m_numIterations = NUM_ITS_INITIAL;
        setDynamic(DYNAMIC_MODE);
        m_viewer = viewer;
    }

    // Changes the behavior of the simulator:
    // If it is set to non-dynamic, it behaves like a ShapeUp style shape defomration
    // framework.
    // If it is set to dynamics, it behaves like a Projective Dynamics simulation.
    void setDynamic(bool dynamic) {
        m_simulator.setDynamic(dynamic);
    }

    // Create the GUI for controlling a Projective Dynamics simulation:
    // starting, pausing and resetting a simulation, as well as adding
    // constraints, defining the number of iterations, etc.
    void initSimulationGUI(bool addConstraintGUI = false) {
        Window* pd_win = new Window(m_viewer, "Simulation Controls");
        pd_win->setPosition(Vector2i(15, 230));
        pd_win->setLayout(new GroupLayout());

        Widget* panel = new Widget(pd_win);
        panel->setLayout(new BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Middle, 0, 10));

        m_startButton = new Button(panel, "Run Simulation");
        m_startButton->setFlags(Button::RadioButton);
        m_startButton->setPushed(false);
        m_startButton->setCallback([this]() {
            start();
        });

        m_stopButton = new Button(panel, "Stop Simulation");
        m_stopButton->setFlags(Button::RadioButton);
        m_stopButton->setPushed(true);
        m_stopButton->setCallback([this]() {
            stop();
        });

        CheckBox* updateNormalsCB = new CheckBox(pd_win, "Update Normals");
        updateNormalsCB->setChecked(m_updateNormals);
        updateNormalsCB->setCallback([this](bool state) {
            m_updateNormals = state;
        });

        Button* reset_b = new Button(pd_win, "Reset Positions");
        reset_b->setCallback([this]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.resetPositions();
            if (was_active) {
                start();
            }
            else {
                uploadPositions();
            }
        });


        Button* addtets_b = new Button(pd_win, "Tetrahedralize");
        addtets_b->setFlags(Button::RadioButton);
        addtets_b->setPushed(false);
        addtets_b->setCallback([this]() {
            is_tetra = true;
            setMesh(true);
        });


        /** NEW! Add UI for temperature diffusion !*/
        PopupButton* popupBtnDiff = new PopupButton(pd_win, "Temp. Diffuse", ENTYPO_ICON_LINK);
        Popup* popupDiff = popupBtnDiff->popup();
        popupDiff->setLayout(new GroupLayout());

        Button* no_diff = new Button(popupDiff, "None");
        Button* uniform = new Button(popupDiff, "Uniform");
        Button* cotan = new Button(popupDiff, "Cotan");

        no_diff->setFlags(Button::RadioButton);
        no_diff->setPushed(true);
        uniform->setFlags(Button::RadioButton);
        cotan->setFlags(Button::RadioButton);

        uniform->setCallback([this, popupBtnDiff, uniform, cotan]() {
            this->diffusion_type = UNILAPLACE;
            diff_activated = uniform->pushed();
            uniform->setPushed(diff_activated);
            cout << "Diffusion Activated: Uniform" << endl;
            popupBtnDiff->setPushed(false);
        });

        cotan->setCallback([this, popupBtnDiff, uniform, cotan]() {
            this->diffusion_type = COTANLAPLACE;
            diff_activated = cotan->pushed();
            cout << "Diffusion Activated: Weighted Cotan" << endl;
            popupBtnDiff->setPushed(false);
        });


        no_diff->setCallback([this, popupBtnDiff, uniform, cotan]() {
            diff_activated = false;
            cout << "Diffusion DeActivated " << endl;
            popupBtnDiff->setPushed(false);
        });

        /** NEW!  change decay **/
        Label* decay_label = new Label(popupDiff, "Decay: ");
        FloatBox<float>* decay_box = new FloatBox<float>(popupDiff);
        decay_box->setEditable(true);
        decay_box->setDefaultValue("0.0");
        decay_box->setMinValue(0.0f);
        decay_box->setMaxValue(1.0f);
        decay_box->setCallback([this](float decay) {
            // Clamp decay between 0.0 and 1.0
            this->decay = max(0.0f, min(1.0f, decay));
        });

        /** NEW!  change number of iterations */
        Label* diff_iterations_label = new Label(popupDiff, "Iterations: ");
        IntBox<int>* diff_iterations_box = new IntBox<int>(popupDiff);
        diff_iterations_box->setEditable(true);
        diff_iterations_box->setDefaultValue("1");
        diff_iterations_box->setMinValue(1);
        diff_iterations_box->setCallback([this](int iterations) {
            this->diffusion_iterations = iterations < 1 ? 1 : iterations;
        });

        /////////////////////////////////////////////////////////////

        PopupButton* popupBtn = new PopupButton(pd_win, "Add constraints", ENTYPO_ICON_LINK);
        Popup* popup = popupBtn->popup();
        popup->setLayout(new GroupLayout());

        Button* b = new Button(popup, "Floor");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            m_simulator.addFloorConstraints(10., 3.);
            m_viewer->setFloorHeight(m_simulator.getFloorHeight());
            m_viewer->showFloor(true);
            updateConstraintsGUI();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Edge Springs");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addEdgeSpringConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Triangle Strain");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addTriangleStrainConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Triangle Bending");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addBendingConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        b = new Button(popup, "Tet Strain");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addTetStrainConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        /** NEW: button for adding temperature elasticity constraints **/
        b = new Button(popup, "Temperature Elasticity");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addEdgeTemperatureElasticityConstraints();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        /** NEW: button for adding position constraints **/
        b = new Button(popup, "Positions");
        b->setCallback([this, popupBtn]() {
            bool was_active = m_simActive;
            stop();
            addPositionConstraintGroup();
            if (was_active) {
                start();
            }
            popupBtn->setPushed(false);
        });

        Button* clear_b = new Button(pd_win, "Clear constraints");
        clear_b->setCallback([this]() {
            stop();
            clearConstraints();
            m_simulator.resetPositions();
            m_viewer->showFloor(false);
            uploadPositions();
            updateConstraintsGUI();
        });

        Label* iterations_label = new Label(pd_win, "Num Loc-Glob Its: ");
        IntBox<int>* iterations_box = new IntBox<int>(pd_win, m_numIterations);
        iterations_box->setEditable(true);
        iterations_box->setCallback([this](int num_its) {
            m_numIterations = num_its;
        });

        if (addConstraintGUI) {
            initConstraintsGUI();
        }

        m_viewer->performLayout();
    }

    // Add an (empty) GUI window in which each constraint group that is added through
    // this API gets a slider which controls a weight multiplier.
    void initConstraintsGUI() {
        m_constraint_window = new Window(m_viewer, "Constraints");
        m_constraint_window->setPosition(Vector2i(700, 25));
        m_constraint_window->setLayout(new GroupLayout());
    }

    // Each time the constraints of the simulation change, this gets called
    // to create one slider (and textboxes) for each constraint group.
    void updateConstraintsGUI() {
        if (!m_constraint_window) return;

        // Clear all previous constraint controls
        while (m_constraint_window->children().size() > 0) {
            m_constraint_window->removeChild(0);
        }

        // For each constraint group add controls
        for (const auto& elem : m_simulator.getConstraintGroups()) {
            const auto g = elem.second;
            new Label(m_constraint_window, g->name, "sans-bold");
            Widget* panel = new Widget(m_constraint_window);
            panel->setLayout(new BoxLayout(nanogui::Orientation::Horizontal, nanogui::Alignment::Middle, 0, 10));

            // Add a sliderand set defaults
            ConstraintSlider* slider = new ConstraintSlider(panel, m_viewer, m_simulator.getNumVerts(), g);

            // Re-initialize system and update positions once the user lets go of the slider
            slider->setFinalCallback([this](float v) {
                bool wasRunning = m_simActive;
                stop();
                m_simulator.initializeSystem();
                update();
                if (wasRunning) start();
            });
        }

        Button* b = new Button(m_constraint_window, "Update");
        b->setCallback([this]() {
            bool wasRunning = m_simActive;
            stop();
            update();
            if (wasRunning) start();
        });

        m_viewer->performLayout();
    }

    // Called when a new mesh is set in the viewer.
    // We convert the Surface_mesh vertices and triangles into Eigen
    // matrices and pass them to the simulator.
    // Additionally, if desired, we generate tetrahedrons that fill the
    // mesh.
    bool setMesh(bool add_tets) {
        std::cout << "New mesh was loaded, re-initializing simulation..." << std::endl;

        // Stop the running simulation
        stop();

        // Convert surface mesh to Eigen matrices
        surface_mesh::Surface_mesh* mesh = m_viewer->getMesh();
        int j = 0;
        ProjDyn::Positions vertices(mesh->n_vertices(), 3);
        ProjDyn::Triangles faces(mesh->n_faces(), 3);
        for (auto f : mesh->faces()) {
            int k = 0;
            for (auto v : mesh->vertices(f)) {
                faces(j, k) = (ProjDyn::Index)v.idx();
                ++k;
            }
            ++j;
        }
        j = 0;
        for (auto v : mesh->vertices()) {
            vertices.row(j) << (ProjDyn::Scalar)mesh->position(v).x,
                (ProjDyn::Scalar)mesh->position(v).y,
                (ProjDyn::Scalar)mesh->position(v).z;
            ++j;
        }

        ProjDyn::Tetrahedrons tets(0, 4);
        if (add_tets) {
            ProjDyn::Positions vol_verts(0, 3);
            try {
                ProjDyn::tetrahedralizeMesh(vertices, faces, vol_verts, tets, 2);
            }
            catch (int e) {
                std::cout << "Error while generating tet-mesh; error-code: " << e << ". Read the console log above for details." << std::endl;
                return false;
            }
            vertices = vol_verts;
            init_temperatures(vertices.size());
        }

        // Set the mesh in the simulator
        m_simulator.setMesh(vertices, faces, tets);

        // Compute neighbourhood info
        m_vertexStars = ProjDyn::makeVertexStars(vertices.size(), faces);

        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(false);
            m_stopButton->setPushed(true);
        }

        updateConstraintsGUI();
        m_viewer->showFloor(false);

        uploadPositions(true);

        return true;
    }

    // Performs a time step and updates the positions that are drawn in the shader window
    bool update(bool forcedUpload = false) {
        if (!m_simulator.isInitialized()) {
            if (!m_simulator.initializeSystem())
                return false;
        }

        /** NEW: update elasticity constraints based on temperature **/
        // 1. Update temperature edge spring and triangle bending constraints
        updateTemperatureConstraints();
        // 2. Recompute LHS, RHS and update solver
        // TODO verify we recompute all required things for triangle bending too
        m_simulator.recomputeConstraintsPrecomputations();
        /** NEW! update temperature values if diffusion activated */
        if(diff_activated) {
            cout << "Diffusing temperatures ";
            switch(diffusion_type) {
                case UNILAPLACE:
                    cout << "using uniform laplacian" << endl;
                    uniform_diffuse();
                    break;
                case COTANLAPLACE:
                    cout << "using cotan weighted laplacian" << endl;
                    weighted_diffuse();
                    break;
                default:
                    break;
            }
        }

        // Simulate one time step
        m_simulator.step(m_numIterations);

        return uploadPositions(forcedUpload);
    }

    // Starts a thread that runs the simulation by constantly
    // calling update(), pausing to not run faster than the set FPS
    bool start() {
        stop();

        // Make sure the simulator is properly initialized
        if (!m_simulator.isInitialized()) {
            if (!m_simulator.initializeSystem()) {
                if (m_startButton && m_stopButton) {
                    m_startButton->setPushed(false);
                    m_stopButton->setPushed(true);
                }
                return false;
            }
        }

        // Create a thread that runs the simulation
        // It calls a function that triggers a time-step every 1000/FPS milliseconds
        m_simActive = true;
        m_simulationThread = std::thread(
            [this]() {
            std::chrono::milliseconds time(1000 / FPS);
            while (m_simActive) {
                std::this_thread::sleep_for(time);
                update();
                glfwPostEmptyEvent();
            }
        }
        );

        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(true);
            m_stopButton->setPushed(false);
        }

        return true;
    }

    // Pauses/stops the current simulation by killing the active thread.
    void stop() {
        if (m_simActive) {
            m_simActive = false;
            m_simulationThread.join();
        }
        if (m_startButton && m_stopButton) {
            m_startButton->setPushed(false);
            m_stopButton->setPushed(true);
        }
    }

    // Extract positions, convert them to column-wise tripples of floats and
    // upload them to the OpenGL buffer
    bool uploadPositions(bool forcedUpload = false) {
        const ProjDyn::Positions& pos = m_simulator.getPositions();

        // Initialize matrix if not done already
        if (m_uploadPos.cols() != m_simulator.getNumOuterVerts() || m_uploadPos.rows() != 3) {
            m_uploadPos.resize(3, m_simulator.getNumOuterVerts());
        }

        // There might be more vertices in the simulation, since inner vertices are not present in the m_viewer
#pragma omp parallel for
        for (int i = 0; i < m_uploadPos.cols(); i++) {
            m_uploadPos(0, i) = (float)pos(i, 0);
            m_uploadPos(1, i) = (float)pos(i, 1);
            m_uploadPos(2, i) = (float)pos(i, 2);
        }
        m_viewer->updateShaderVertices(m_uploadPos, forcedUpload);

        if (m_updateNormals && m_vertexStars.size() >= m_simulator.getNumOuterVerts()) {
            // Initialize matrix if not done already
            if (m_uploadNormals.cols() != pos.rows() || m_uploadNormals.rows() != 3) {
                m_uploadNormals.resize(3, m_simulator.getNumOuterVerts());
            }
            // Compute per-triangle normals and sum them on vertices
            const ProjDyn::Triangles& tris = m_simulator.getTriangles();
            m_uploadNormals.setZero();
#pragma omp parallel for
            for (int vInd = 0; vInd < m_simulator.getNumOuterVerts(); vInd++) {
                float fac = 1.f / m_vertexStars[vInd].size();
                ProjDyn::Vector3 normal;
                normal.setZero();
                for (int locInd = 0; locInd < m_vertexStars[vInd].size(); locInd++) {
                    int t = m_vertexStars[vInd][locInd].t1;
                    normal += fac * (pos.row(tris(t, 2)) - pos.row(tris(t, 1))).cross((pos.row(tris(t, 0)) - pos.row(tris(t, 1)))).normalized();
                }
                m_uploadNormals(0, vInd) = normal(0);
                m_uploadNormals(1, vInd) = normal(1);
                m_uploadNormals(2, vInd) = normal(2);
            }

            //Upload the normals
            m_viewer->updateShaderNormals(m_uploadNormals, forcedUpload);
        }

        return true;
    }

    // Update elasticity constraints based on temperature
    void updateTemperatureConstraints() {
        Surface_mesh* mesh = m_viewer->getMesh();
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        Surface_mesh::Vertex_property<Scalar> v_temperature = mesh->vertex_property<Scalar>("v:temperature", 0.0);
        const auto v_lookup_table = m_viewer->getVertexLookupTable();
        const auto groups = m_simulator.getConstraintGroups();

        // Find group of temperature based constraints, return false if it doesn't exist
        auto elem = groups.find("Edge Temperature Elasticity");
        if (elem != groups.end()) {
            auto cg = elem->second;
            // Recompute weight of each constraint of the group and update
            for (const auto c : cg->constraints) {
                const std::vector<Index>& vIndices = c->getIndices();

                // Each vertex may either be on the surface or in the interior
                const auto t0 = vIndices[0] < mesh->n_vertices() ?
                    v_temperature[v_lookup_table[vIndices[0]]] : m_temperatures[vIndices[0]];
                const auto t1 = vIndices[1] < mesh->n_vertices() ?
                    v_temperature[v_lookup_table[vIndices[1]]] : m_temperatures[vIndices[1]];

                const Scalar avgTemp = 0.5 * (t0 + t1);
                if (avgTemp >= 200) {
                    c->setWeight(0.0001f);
                } else {
                    const Scalar edgeLen = (sim_verts.row(vIndices[0]) - sim_verts.row(vIndices[1])).norm();
                    c->setWeight(edgeLen);
                }
            }
        } else {
            std::cout << "Warning: no temperature elasticity constraint group found." << std::endl;
        }

        // Find group of triangle bending constraints
        elem = groups.find("Tri Bending");
        if (elem != groups.end()) {
            const auto cg = elem->second;
            // Recompute weight of each constraint in the group
            for (const auto c: cg->constraints) {
                const std::vector<Index>& vIndices = c->getIndices();
                const Scalar edgeLen = (sim_verts.row(vIndices[0]) - sim_verts.row(vIndices[1])).norm();
                const auto v0 = v_lookup_table[vIndices[0]];
                const auto v1 = v_lookup_table[vIndices[1]];

                const Scalar t0 = v_temperature[v0];
                const Scalar t1 = v_temperature[v1];
                const Scalar avgTemp = 0.5 * (t0 + t1);

                if (avgTemp >= 200) {
                    c->setWeight(0.001f);
                } else {
                    c->setWeight(edgeLen);
                }
            }
        } else {
            std::cout << "Warning: no temperature elasticity constraint group found." << std::endl;
        }
    }

    // Set external forces to point into downwards y direction with a certain magnitude
    void setGravity(ProjDyn::Scalar g) {
        m_simulator.setGravity(g);
    }

    // Can be used to visually mark certain vertices based on an integer
    // value they receive in a row-vector (vertex ID corresponds to column)
    void uploadVertexStatus(const Eigen::Matrix<int, 1, -1>& vStatus) {
        m_viewer->updateVertexStatus(vStatus);
    }

    // Change the number of local-global iterations in the simulator
    void setNumIterations(Index numIts) {
        m_numIterations = numIts;
    }
    Index getNumIterations() {
        return m_numIterations;
    }

    // Called when a new mesh was set in the viewer
    bool setMesh() {
        return setMesh(false);
    }

    // Add constraints to the simulator, either as groups, as lists
    // or as single constraints 
    void addConstraints(const std::vector<ProjDyn::ConstraintPtr>& constraints) {
        m_simulator.addConstraints(constraints);
        updateConstraintsGUI();
    }
    void addConstraints(ProjDyn::ConstraintGroupPtr constraints) {
        m_simulator.addConstraints(constraints);
        updateConstraintsGUI();
    }
    void addConstraint(const ProjDyn::ConstraintPtr& constraint) {
        m_simulator.addConstraint(constraint);
        updateConstraintsGUI();
    }

    // Remove a constraint group from the simulator, either by value or by name.
    void removeConstraints(ProjDyn::ConstraintGroupPtr constraints) {
        m_simulator.removeConstraints(constraints);
        updateConstraintsGUI();
    }
    void removeConstraints(std::string constraint_name) {
        m_simulator.removeConstraints(constraint_name);
        updateConstraintsGUI();
    }

    // Remove all constraints from the simulation.
    void clearConstraints() {
        m_simulator.clearConstraints();
    }

    // In the following are several helper function to add several types of specific
    // constraints to the simulation:


    // Add tetrahedral strain constraints to all tets:
    void addTetStrainConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allTets;
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        for (Index i = 0; i < tets.rows(); i++) allTets.push_back(i);
        addTetStrainConstraints(allTets, weight);
    }

    // Add tetrahedral strain constraints to some tets:
    void addTetStrainConstraints(const std::vector<Index>& tetInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> tet_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        if (tets.rows() > 0) {
            for (Index i : tetInds) {
                if (i >= tets.rows()) continue;
                std::vector<ProjDyn::Index> tetInds;
                for (int j = 0; j < 4; j++) tetInds.push_back(tets(i, j));
                // The weight is the tet volume
                ProjDyn::Scalar w = ProjDyn::tetrahedronVolume(sim_verts, tets.row(i));
                if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                    ProjDyn::TetStrainConstraint* esc = new ProjDyn::TetStrainConstraint(tetInds, w, sim_verts);
                    tet_constraints.push_back(std::shared_ptr<ProjDyn::TetStrainConstraint>(esc));
                }
            }
            addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tet Strain", tet_constraints, weight));
        }
        else {
            std::cout << "WARNING: Could not add tet-strain - no tets available!" << std::endl;
        }
    }

    // Add bending constraints to all vertices of the simulation
    void addBendingConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allOuterVerts;
        for (Index i = 0; i < m_simulator.getNumOuterVerts(); i++) allOuterVerts.push_back(i);
        addBendingConstraints(allOuterVerts, weight);
    }

    // Add bending constraints to some vertices of the simulation
    void addBendingConstraints(const std::vector<Index>& vertInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> bend_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        std::vector<ProjDyn::VertexStar>& vStars = m_vertexStars;
        ProjDyn::Vector voronoiAreas = ProjDyn::vertexMasses(m_simulator.getPositions(), m_simulator.getTriangles());
        Surface_mesh* smesh = m_viewer->getMesh();
        for (auto v : smesh->vertices()) {
            Index i = v.idx();
            if (std::find(vertInds.begin(), vertInds.end(), i) == vertInds.end()) continue;
            if (smesh->is_boundary(v)) continue;
            if (i >= m_simulator.getNumOuterVerts()) continue;
            // The weight is the voronoi area
            ProjDyn::Scalar w = voronoiAreas(i) * 0.01;
            if (w > 1e-6) {
                // The constraint is constructed, made into a shared pointer and appended to the list
                // of constraints.
                ProjDyn::BendingConstraint* esc = new ProjDyn::BendingConstraint(vStars[i], w, voronoiAreas(i), sim_verts, m_simulator.getTriangles());
                bend_constraints.push_back(std::shared_ptr<ProjDyn::BendingConstraint>(esc));
            }
        }
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tri Bending", bend_constraints, weight));
    }

    // Add triangular strain constraints to all triangles of the simulation
    void addTriangleStrainConstraints(ProjDyn::Scalar weight = 1.) {
        std::vector<Index> allTris;
        for (Index i = 0; i < m_simulator.getTriangles().rows(); i++) allTris.push_back(i);
        addTriangleStrainConstraints(allTris, weight);
    }

    // Add triangular strain constraints to some triangles of the simulation
    void addTriangleStrainConstraints(const std::vector<Index>& triInds, ProjDyn::Scalar weight = 1.) {
        std::vector<ProjDyn::ConstraintPtr> tri_constraints;
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        const ProjDyn::Triangles& tris = m_simulator.getTriangles();
        for (Index i : triInds) {
            if (i >= tris.rows()) continue;
            std::vector<ProjDyn::Index> triInds;
            for (int j = 0; j < 3; j++) triInds.push_back(tris(i, j));
            // The weight is the triangle area
            ProjDyn::Scalar w = ProjDyn::triangleArea(sim_verts, tris.row(i));
            if (w > 1e-6) {
                // The constraint is constructed, made into a shared pointer and appended to the list
                // of constraints.
                ProjDyn::TriangleStrainConstraint* esc = new ProjDyn::TriangleStrainConstraint(triInds, w, sim_verts);
                tri_constraints.push_back(std::shared_ptr<ProjDyn::TriangleStrainConstraint>(esc));
            }
        }
        // Add constraints
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Tri Strain", tri_constraints, weight));
    }

    // Add edge spring constraints to all edges of the simulation
    void addEdgeSpringConstraints(ProjDyn::Scalar weight = 1.) {
        // For tet meshes we cannot use Surface_mesh
        if (m_simulator.getTetrahedrons().rows() > 0) {
            addEdgeSpringConstraintsTets(weight);
        }
        else {
            const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
            const ProjDyn::Triangles& tris = m_simulator.getTriangles();
            Surface_mesh* smesh = m_viewer->getMesh();
            std::vector<ProjDyn::ConstraintPtr> spring_constraints;
            for (auto edge : smesh->edges()) {
                // The weight is set to the edge length
                ProjDyn::Scalar w = (sim_verts.row(smesh->vertex(edge, 0).idx()) - sim_verts.row(smesh->vertex(edge, 1).idx())).norm();
                if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                    std::vector<Index> edge_inds;
                    edge_inds.push_back(smesh->vertex(edge, 0).idx());
                    edge_inds.push_back(smesh->vertex(edge, 1).idx());
                    ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge_inds, w, sim_verts);
                    spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                }
            }
            addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Springs", spring_constraints, weight));
        }
    }

    /** NEW: Add edge springs based on temperature of each edge (avg temp of both vertices) **/
    void addEdgeTemperatureElasticityConstraints(ProjDyn::Scalar weight = 1.) {
        if (m_simulator.getTetrahedrons().rows() > 0) {
            addEdgeTemperatureSpringConstraintsTets(weight);
        }
        else {
            const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
            const ProjDyn::Triangles& tris = m_simulator.getTriangles();
            Surface_mesh* smesh = m_viewer->getMesh();
            Surface_mesh::Vertex_property<Scalar> v_temperature = smesh->vertex_property<Scalar>("v:temperature", 0.0);
            std::vector<ProjDyn::ConstraintPtr> spring_constraints;
            for (auto edge : smesh->edges()) {
                /** NEW: w = edge_len / (1 + avg.temp) **/
                ProjDyn::Scalar w = (sim_verts.row(smesh->vertex(edge, 0).idx()) - sim_verts.row(smesh->vertex(edge, 1).idx())).norm();
                const Scalar t0 = v_temperature[smesh->vertex(edge, 0)];
                const Scalar t1 = v_temperature[smesh->vertex(edge, 1)];
                const Scalar edgeTemp = 0.5 * (t0 + t1);
                w /= (1.0 + edgeTemp);

                if (w > 1e-6) {
                    // The constraint is constructed, made into a shared pointer and appended to the list
                    // of constraints.
                    std::vector<Index> edge_inds;
                    edge_inds.push_back(smesh->vertex(edge, 0).idx());
                    edge_inds.push_back(smesh->vertex(edge, 1).idx());
                    ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge_inds, w, sim_verts);
                    spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                }
            }
            addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Temperature Elasticity", spring_constraints, weight));
        }
    }

    void addPositionConstraintGroup(ProjDyn::Scalar weight= 1.0f) {
        auto selected_verts = m_viewer->getSelectedVertices();
        if(selected_verts.empty()) {
            cout << "Please select vertices to fix positions" << endl;
            return;
        }
        const ProjDyn::Positions& curPos = getPositions();
        ProjDyn::Positions con_pos;
        con_pos.setZero(selected_verts.size(), 3);
        Index ind = 0;
        for (Index v : selected_verts) {
            con_pos.row(ind) = curPos.row(v);
            ind++;
        }
        // Create constraint groups and add them to the simulation
        std::shared_ptr<ProjDyn::PositionConstraintGroup> con = std::shared_ptr<ProjDyn::PositionConstraintGroup>(new ProjDyn::PositionConstraintGroup(selected_verts, weight, con_pos));

        auto conGroup = std::make_shared<ProjDyn::ConstraintGroup>("Fixed Pos.", std::vector<ProjDyn::ConstraintPtr>({ con }), 1.);
        addConstraints(conGroup);
    }

    // Gets called when the user is grabbing some vertices with the mouse
    // (alt + left-click).
    void grab(const std::vector<ProjDyn::Index>& grabbedVerts, const std::vector<Vector3f>& grabPos) {
        if (!m_simulator.isInitialized() || !m_simActive) return;

        m_simulator.setGrab(grabbedVerts, grabPos);
    }

    // Gets called when the user is letting go of the grabbed vertices.
    void releaseGrab() {
        if (!m_simulator.isInitialized()) return;

        m_simulator.releaseGrab();
    }

    // In the following are some simple getter/setter function to communicate
    // with the simulation:

    const ProjDyn::Positions& getPositions() {
        return m_simulator.getPositions();
    }

    const ProjDyn::Triangles& getTriangles() {
        return m_simulator.getTriangles();
    }

    const ProjDyn::Tetrahedron& getTetrahedra() {
        return m_simulator.getTetrahedrons();
    }

    const ProjDyn::Index getNumVertices() {
        return m_simulator.getNumVerts();
    }

private:
    // Global variables used by the callback functions
    int m_numIterations = NUM_ITS_INITIAL;
    ProjDyn::Simulator m_simulator;
    Eigen::MatrixXf m_uploadPos, m_uploadNormals;
    std::thread m_simulationThread;
    bool m_simActive = false;
    std::vector<ProjDyn::VertexStar> m_vertexStars;
    Window* m_constraint_window = nullptr;
    Viewer* m_viewer;
    Button* m_startButton = nullptr;
    Button* m_stopButton = nullptr;
    bool m_updateNormals = UPDATE_NORMALS;

    void addEdgeSpringConstraintsTets(ProjDyn::Scalar weight = 1.) {
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        std::vector<ProjDyn::ConstraintPtr> spring_constraints;
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        // If tets are available, add a spring on each tet-edge
        for (Index i = 0; i < tets.rows(); i++) {
            for (int j = 0; j < 4; j++) {
                std::vector<ProjDyn::Index> edge;
                edge.push_back(tets(i, j));
                edge.push_back(tets(i, (j + 1) % 4));
                // Easy way to make sure each edge only gets added once:
                if (edge[0] < edge[1]) {
                    // The weight is set to the edge length
                    ProjDyn::Scalar w = (sim_verts.row(edge[0]) - sim_verts.row(edge[1])).norm();
                    if (w > 1e-6) {
                        // The constraint is constructed, made into a shared pointer and appended to the list
                        // of constraints.
                        ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge, w, sim_verts);
                        spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                    }
                }
            }
        }

        // Add constraints
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Springs", spring_constraints, weight));
    }

    void addEdgeTemperatureSpringConstraintsTets(ProjDyn::Scalar weight = 1.) {
        const ProjDyn::Positions& sim_verts = m_simulator.getInitialPositions();
        std::vector<ProjDyn::ConstraintPtr> spring_constraints;
        const ProjDyn::Tetrahedrons& tets = m_simulator.getTetrahedrons();
        // If tets are available, add a spring on each tet-edge
        for (Index i = 0; i < tets.rows(); i++) {
            for (int j = 0; j < 4; j++) {
                std::vector<ProjDyn::Index> edge;
                edge.push_back(tets(i, j));
                edge.push_back(tets(i, (j + 1) % 4));
                // Easy way to make sure each edge only gets added once:
                if (edge[0] < edge[1]) {
                    // The weight is set to the edge length
                    ProjDyn::Scalar w = (sim_verts.row(edge[0]) - sim_verts.row(edge[1])).norm();
                    if (w > 1e-6) {
                        // The constraint is constructed, made into a shared pointer and appended to the list
                        // of constraints.
                        ProjDyn::EdgeSpringConstraint* esc = new ProjDyn::EdgeSpringConstraint(edge, w, sim_verts);
                        spring_constraints.push_back(std::shared_ptr<ProjDyn::EdgeSpringConstraint>(esc));
                    }
                }
            }
        }

        // Add constraints
        addConstraints(std::make_shared<ProjDyn::ConstraintGroup>("Edge Temperature Elasticity", spring_constraints, weight));
    }

    /** Start of Mesh Processing functions
    //todo create class to harbor all processing

     *  Weight calculations are functions taken from:
     *  //=============================================================================
        //
        //   Code framework for the lecture
        //
        //   "Digital 3D Geometry Processing"
        //
        //   Gaspard Zoss, Alexandru Ichim
        //
        //   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
        //         EPF Lausanne
        //
        //-----------------------------------------------------------------------------

     *  */
    void calc_weights() {
        calc_edges_weights();
        calc_vertices_weights();
    }

    void calc_edges_weights() {
        auto mesh_ = m_viewer->getMesh();
        auto e_weight = mesh_->edge_property<Scalar>("e:weight", 0.0f);
        auto points = mesh_->vertex_property<Point>("v:point");

        Surface_mesh::Halfedge h0, h1, h2;
        Point p0, p1, p2, d0, d1;

        for (auto e: mesh_->edges())
        {
            e_weight[e] = 0.0;

            h0 = mesh_->halfedge(e, 0);
            p0 = points[mesh_->to_vertex(h0)];

            h1 = mesh_->halfedge(e, 1);
            p1 = points[mesh_->to_vertex(h1)];

            if (!mesh_->is_boundary(h0))
            {
                h2 = mesh_->next_halfedge(h0);
                p2 = points[mesh_->to_vertex(h2)];
                d0 = p0 - p2;
                d1 = p1 - p2;
                e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
            }

            if (!mesh_->is_boundary(h1))
            {
                h2 = mesh_->next_halfedge(h1);
                p2 = points[mesh_->to_vertex(h2)];
                d0 = p0 - p2;
                d1 = p1 - p2;
                e_weight[e] += dot(d0,d1) / norm(cross(d0,d1));
            }
        }
    }

    void calc_vertices_weights() {
        Surface_mesh::Face_around_vertex_circulator vf_c, vf_end;
        Surface_mesh::Vertex_around_face_circulator fv_c;
        Scalar area;
        auto mesh_ = m_viewer->getMesh();
        auto v_weight = mesh_->vertex_property<Scalar>("v:weight", 0.0f);

        for (auto v: mesh_->vertices()) {
            area = 0.0;
            vf_c = mesh_->faces(v);

            if(!vf_c) {
                continue;
            }

            vf_end = vf_c;

            do {
                fv_c = mesh_->vertices(*vf_c);

                const Point& P = mesh_->position(*fv_c);  ++fv_c;
                const Point& Q = mesh_->position(*fv_c);  ++fv_c;
                const Point& R = mesh_->position(*fv_c);

                area += norm(cross(Q-P, R-P)) * 0.5f * 0.3333f;

            } while(++vf_c != vf_end);

            v_weight[v] = 0.5 / area;
        }
    }

    /**
     * Enum type defining the different possible laplacian types
     */
    enum struct LaplacianType{
        CONSTANT=1,
        COTAN=2
    };

    /**
     * Computes the update of ti, following the laplazian expression ti = ti + sum (tj-ti)*wi
     * @param target The target vertex i, for which the temperature ti is to be updated
     * @param mesh The mesh of target
     * @param type Laplacian type. Influences the weighting.
     * @param e_weight Edge weights. Useful in case of cotan Laplacian, ignored for constant (uniform) Laplace.
     * @param temperature Temperature array, per vertex
     * @return New value of the scalar ti for the target vertex
     */
    double normalizedTemperatureLaplacian(const Surface_mesh::Vertex &target, const Surface_mesh &mesh, const LaplacianType type,
                                          const Surface_mesh::Edge_property<Scalar> &e_weight, const Surface_mesh::Vertex_property<Scalar> &temperature){

        Surface_mesh::Halfedge_around_vertex_circulator circulator(&mesh, target);
        double total(0);
        double weight_sum(0), weight(0), diff_timestep(0.499999);
        double target_temperature = temperature[target];
        if(circulator){
            for(auto half_edge: circulator){
                switch(type){
                    case LaplacianType ::CONSTANT:
                        weight = 1.0;
                        break;
                    case LaplacianType ::COTAN:
                        weight=e_weight[mesh.edge(half_edge)];
                        break;
                }
                weight_sum += weight;
                total += (temperature[mesh.to_vertex(half_edge)] - target_temperature) * weight;
            }
        }

        if(weight_sum != 0.0){
            total /= weight_sum;
        }

        return target_temperature + diff_timestep * total;
    }

    /**
     * Updates temperatures of the provided mesh, according to a specific Laplacian rule update.
     * @param mesh The target mesh
     * @param v_new_temp Temporary array to hold the temperatures before update
     * @param type Laplacian update type
     * @param e_weight Edge weights of the mesh
     * @param is_source Array property specifying which vertex of the mesh is a source
     * @param temperatures Temperatures of each vertex
     */
    void applySmoothing(const Surface_mesh &mesh, Surface_mesh::Vertex_property<Scalar> &v_new_temp, const LaplacianType &type, const Surface_mesh::Edge_property<Scalar> &e_weight,
            const Surface_mesh::Vertex_property<bool> &is_source, Surface_mesh::Vertex_property<Scalar> &temperatures){
        for(auto vertex : mesh.vertices()) {
            if(!mesh.is_boundary(vertex) && !is_source[vertex]){
                switch(type){
                    case LaplacianType::CONSTANT:
                        v_new_temp[vertex]= normalizedTemperatureLaplacian(vertex, mesh, LaplacianType::CONSTANT,e_weight, temperatures);
                        break;
                    case LaplacianType ::COTAN:
                        v_new_temp[vertex]= normalizedTemperatureLaplacian(vertex, mesh, LaplacianType::COTAN,e_weight, temperatures);
                        break;
                }
            }
        }

        for(auto vertex: mesh.vertices()){
            if(!mesh.is_boundary(vertex) && !is_source[vertex]){
                temperatures[vertex] = v_new_temp[vertex];
            }
        }
    }

    /**
     * If the m_temperatures array was not initialized yet, initialize it.
     * The way this is done is by creating a vector of size n_vertices. The first n values (n the number of vertices on the surface of the mesh)
     *  are copied back from the v:temperature vertex property. Other values are simply initialized to 0.
     *  If the mesh was initialized, the temperatures of the envelope are still copied back (in case user updated them). Others are unchanged.
     * @param n_vertices Number of vertices in the tetrahedralized mesh.
     */
    void init_temperatures(const unsigned int n_vertices){
        auto mesh_ = m_viewer->getMesh();
        Eigen::VectorXd tmp_temps(n_vertices);
        Surface_mesh::Vertex_property<Scalar> v_temp = mesh_->vertex_property<Scalar>("v:temperature",0.0);

        // Store fist these temperatures
        Index i = 0;
        for(auto v: mesh_->vertices()){
            i = v.idx();
            tmp_temps(i) = v_temp[v];
        }
        if(!init_temp){

            for(Index j = i+1; j < n_vertices; ++j){
                tmp_temps(j) = 0.0;
            }

            init_temp = true;
        } else {
            for(Index j = i+1; j < n_vertices; ++j){
                tmp_temps(j) = m_temperatures(j);
            }
        }

        m_temperatures = tmp_temps;
    }

    /**
     * Performs cotan weighted Laplacian of the temperatures
     * Updates the mesh's temperatures according to the heat equation, using cotangent weights.
     * The cotan weighted Laplacian is either surface-based for regular meshes, or volumetric, in case of tetrahedralized meshes.
     * Once new temperatures have been computed, it also updates the displayed temperatures for visualization.
     */
    void weighted_diffuse() {
        cout << "Weighted diffusion with decay value : " << decay << " for " << diffusion_iterations << " iterations" <<  endl;

        Surface_mesh::Vertex_around_vertex_circulator vv_c, vv_end;
        double laplacian;
        unsigned int w;
        auto mesh_ = m_viewer->getMesh();

        const double one_ov_four(1.0/4.0);
        const double one_ov_six(1.0/6.0);
        if(is_tetra){
            /// This code solves the equation (D^-1  - delta * L) P(t+1) = D^-1 P(t),
            /// for P(t+1). L is the matrix of cotangent weights. D^-1 is the mass-matrix of volumes (denoted M from here on).
            /// P(t) is the matrix of temperatures.

            const ProjDyn::Tetrahedrons &tets = m_simulator.getTetrahedrons();
            const ProjDyn::Positions &pos = m_simulator.getPositions();

            /// Computes the cotangent matrix L
            Eigen::SparseMatrix<double> L;
            igl::cotmatrix(pos, m_simulator.getTetrahedrons(), L);

            unsigned int n_vertices = L.rows();

            /// If the full temperatures (i.e: for all vertices in the volume) are not initialized, this function will
            /// ensure it is the case, while preserving temperatures of the already existing vertices.
            /// Otherwise it does nothing.
            init_temperatures(n_vertices);

            /// We must now compute M , which is our only missing quantity.
            /// We will also compute M P(t), while we're at it, and store it as B=M P(t).
            std::vector<Eigen::Triplet<double> > triplets;
            Eigen::MatrixXd B(n_vertices, 1);

            unsigned int n_tetra = tets.rows();

            double volumes[n_vertices];

            // Just making sure volumes is zero-initialized (should be per C++ standard, but you never know)
            for(Index j = 0; j < n_vertices; ++j){
                volumes[j] = 0.0;
            }

            /// Fills in the volumes array, which sums up volume contributions of all tetrahedra containing a vertex i.
            /// The volume of a tetrahedron is given by abs((a-d).((b-d) x (c-d)))/6
            for(Index j = 0; j < n_tetra; ++j){
                /// First get index of all vertices of tetrahedron
                Index ind_a = tets(j,0);
                Index ind_b = tets(j,1);
                Index ind_c = tets(j,2);
                Index ind_d = tets(j,3);

                /// Then, extract the positions of all vertices
                Vector<Scalar,3> a(pos(ind_a,0),
                        pos(ind_a,1),
                        pos(ind_a,2));
                Vector<Scalar,3> b(pos(ind_b,0),
                                  pos(ind_b,1),
                                  pos(ind_b,2));
                Vector<Scalar,3> c(pos(ind_c,0),
                                  pos(ind_c,1),
                                  pos(ind_c,2));
                Vector<Scalar,3> d(pos(ind_d,0),
                                  pos(ind_d,1),
                                  pos(ind_d,2));

                /// Finally compute the volume of tetrahedron
                double volume = abs(dot(a-d, cross((b-d),(c-d))))*one_ov_six;

                /// The sum for all vertices must be updated by the current volume
                volumes[ind_a] += volume;
                volumes[ind_b] += volume;
                volumes[ind_c] += volume;
                volumes[ind_d] += volume;
            }

            /// The mass matrix M is computed, as well as the righthand side M*P^(t)=B
            for (Index j = 0; j < n_vertices; ++j) {
                double mii = volumes[j]*one_ov_four;
                triplets.push_back(Eigen::Triplet<double>(j, j, mii));
                double tmp = mii*m_temperatures(j);
                assert(!isnan(tmp)); // NaN values are unacceptable here, as volume should always be defined along with temperature.
                B(j) = tmp;
            }

            /// The mass matrix M should be sparse (more efficient)
            Eigen::SparseMatrix<double> M(n_vertices, n_vertices); // The reason we store in a sparse matrix is to enable operations with L per Eigen specifications
            M.setFromTriplets(triplets.begin(), triplets.end());

            double diff_timestep(0.0001);

            const auto &S = (M - diff_timestep * L); // Because all quantities are sparse, we can directly perform the following !
            Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(S);
            assert(solver.info() == Eigen::Success);

            /// The solver of the form Ax = B received A = (M - delta L). We give it B= M P(t), to retrieve at last P(t+1)
            Eigen::MatrixXd P = solver.solve(B).eval(); // these are the new temperatures!

            /// All that remains is to update the temperatures as well as the displayed temperatures!
            Surface_mesh::Vertex_property<Scalar> v_temp = mesh_->vertex_property<Scalar>("v:temperature", 0.0);
            auto v = mesh_->vertices().begin();
            for (Index i = 0; i < n_vertices; ++i) {
                double r = P(i, 0);
                /// Should temperatures somehow become nan, it means the solver failed, which happens if S is not positive definite.
                /// It was observed that this consistently happens after thousands of iterations where the system is more or less settled.
                /// Therefore, it could be numerical errors piling up until the system is no more defined. For this reason, if a nan is found
                /// it is simply ignored, so that the system remains stable.
                if (!isnan(r)) {
                    m_temperatures(i) = r;
                    if (v != mesh_->vertices().end()) {
                        v_temp[*v] = r;
                        assert(!isnan(v_temp[*v]));
                        ++v;
                    }
                } else {
                    cout << "Reached NaN value (matrix not positive definite)" << endl;
                }
            }
        }else {

            Surface_mesh::Vertex_property<Scalar> v_new_temp = mesh_->vertex_property<Scalar>("v:new_temperatures");
            Surface_mesh::Vertex_property<bool> v_is_source = mesh_->vertex_property<bool>("v:is_source", false);
            Surface_mesh::Vertex_property<Scalar> v_temp = mesh_->vertex_property<Scalar>("v:temperature", 0.0);
            Surface_mesh::Edge_property<Scalar> e_weight = mesh_->edge_property<Scalar>("e:weight", 0.0f);
            for (unsigned int iter = 0; iter < diffusion_iterations; ++iter) {
                // For each non-boundary vertex, update its temperature according to the cotan Laplacian operator
                calc_weights();
                applySmoothing(*mesh_, v_new_temp, LaplacianType::COTAN, e_weight, v_is_source, v_temp);
            }
        }

        m_viewer->setTemperatureColor();

    }


    /**
     * Performs uniform Laplacian update of the temperatures of the mesh
     */
    void uniform_diffuse() {
        Surface_mesh::Vertex_around_vertex_circulator vv_c, vv_end;
        double laplacian;
        unsigned int w;
        auto mesh_ = m_viewer->getMesh();
        double one_ov_four=1.0/4.0;

        if(is_tetra) {
            throw std::invalid_argument("Not implemented for simple case (try cotan instead)");
        } else {
            Surface_mesh::Vertex_property<Scalar> v_new_temp = mesh_->vertex_property<Scalar>("v:new_temperatures");
            Surface_mesh::Vertex_property<bool> v_is_source = mesh_->vertex_property<bool>("v:is_source", false);
            Surface_mesh::Vertex_property<Scalar> v_temp = mesh_->vertex_property<Scalar>("v:temperature",0.0);

            Surface_mesh::Edge_property<Scalar> e_weight;
            for (unsigned int iter=0; iter<diffusion_iterations; ++iter) {
                // For each non-boundary vertex, update its temperature according to the uniform Laplacian operator
                applySmoothing(*mesh_, v_new_temp, LaplacianType::CONSTANT, e_weight, v_is_source, v_temp);
            }
        }

        m_viewer->setTemperatureColor();
    }

};